#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simulador Celeste 3D — somente biblioteca padrão, Python 3.9+.

Executar: python simulador_celeste.py
Testar:   python simulador_celeste.py --testar
Modelo: gravitação newtoniana de N corpos, velocity Verlet em 3 dimensões.
Unidades internas: massa kg, posição UA, velocidade UA/dia, tempo dia.
Tkinter só é importado ao abrir a interface; os testes não exigem tela.
"""

import argparse
from collections import deque
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
import re
import sys
import time
import uuid

AU_KM = 149597870.7
DAY_S = 86400.0
G_SI = 6.67430e-11
G = G_SI * DAY_S**2 / (AU_KM * 1000.0)**3
SUN_KG = 1.98847e30
EARTH_KG = 5.9722e24
KM_S_TO_AU_DAY = DAY_S / AU_KM
KINDS = ("Livre", "Estrela", "Em órbita")
MAX_BODIES = 150


def add(a, b):
    return [a[k] + b[k] for k in range(3)]


def sub(a, b):
    return [a[k] - b[k] for k in range(3)]


def mul(a, x):
    return [v * x for v in a]


def dot(a, b):
    return sum(a[k] * b[k] for k in range(3))


def norm(a):
    return math.sqrt(dot(a, a))


def number(value, name, lower=-1e100, upper=1e100):
    """Aceita vírgula decimal ou ponto, inclusive notação científica."""
    if isinstance(value, bool):
        raise ValueError("%s deve ser um número." % name)
    try:
        result = float(str(value).strip().replace(",", "."))
    except (ValueError, TypeError):
        raise ValueError("%s: informe um número válido." % name) from None
    if not math.isfinite(result) or not lower <= result <= upper:
        raise ValueError("%s deve estar entre %g e %g." % (name, lower, upper))
    return result


@dataclass
class Body:
    name: str
    mass: float
    radius_km: float
    pos: list
    vel: list
    kind: str = "Livre"
    parent: str = ""
    color: str = "#8dccff"
    visual_px: float = 7.0
    distance_factor: float = 1.0
    uid: str = field(default_factory=lambda: uuid.uuid4().hex)

    def to_dict(self):
        return dict(id=self.uid, nome=self.name, tipo=self.kind,
                    estrela_id=self.parent, massa_kg=self.mass,
                    raio_km=self.radius_km, posicao_ua=list(self.pos),
                    velocidade_km_s=mul(self.vel, 1.0 / KM_S_TO_AU_DAY),
                    cor=self.color, raio_visual_px=self.visual_px,
                    fator_distancia_visual=self.distance_factor)

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict):
            raise ValueError("Cada corpo deve ser um objeto JSON.")
        name = data.get("nome", "")
        kind = data.get("tipo", "Livre")
        color = data.get("cor", "#8dccff")
        uid = data.get("id", "")
        parent = data.get("estrela_id", "")
        if not isinstance(name, str) or not name.strip() or len(name) > 60:
            raise ValueError("Nome do corpo: use de 1 a 60 caracteres.")
        if kind not in KINDS:
            raise ValueError("Tipo de corpo desconhecido.")
        if not isinstance(uid, str) or not uid or len(uid) > 100:
            raise ValueError("Identificador de corpo inválido.")
        if not isinstance(parent, str):
            raise ValueError("Referência de estrela inválida.")
        if not isinstance(color, str) or not re.fullmatch(r"#[0-9a-fA-F]{6}", color):
            raise ValueError("Cor inválida. Use #RRGGBB.")
        vectors = []
        for key in ("posicao_ua", "velocidade_km_s"):
            values = data.get(key)
            if not isinstance(values, list) or len(values) != 3:
                raise ValueError("%s precisa de três componentes." % key)
            vectors.append([number(v, key, -1e12, 1e12) for v in values])
        return cls(name.strip(), number(data.get("massa_kg"), "Massa", 1e-6, 1e36),
                   number(data.get("raio_km"), "Raio físico", 0, 1e12),
                   vectors[0], mul(vectors[1], KM_S_TO_AU_DAY), kind, parent,
                   color, number(data.get("raio_visual_px", 7), "Raio visual", 1, 200),
                   number(data.get("fator_distancia_visual", 1), "Fator visual", 1e-6, 1e6), uid)


def validate_links(bodies):
    if len(bodies) > MAX_BODIES:
        raise ValueError("Limite de %d corpos por cenário." % MAX_BODIES)
    by_id = {b.uid: b for b in bodies}
    if len(by_id) != len(bodies):
        raise ValueError("Há identificadores de corpos repetidos.")
    for b in bodies:
        if b.kind == "Em órbita":
            p = by_id.get(b.parent)
            if p is None or p.kind != "Estrela" or p.uid == b.uid:
                raise ValueError("%s precisa de uma estrela válida." % b.name)
        elif b.parent:
            raise ValueError("Somente corpos 'Em órbita' têm estrela associada.")


def orbital_state(star, mass, a, e, inclination, node, periapsis, anomaly):
    """Elementos keplerianos elípticos -> estado instantâneo no referencial global.

    a em UA; ângulos em graus; e entre 0 e 0.99. Inclui a massa do corpo
    no parâmetro gravitacional. Outros corpos perturbarão essa órbita.
    """
    a = number(a, "Semieixo maior", 1e-9, 1e9)
    e = number(e, "Excentricidade", 0, 0.99)
    inc, ascending, arg, nu = [math.radians(number(x, "Ângulo", -36000, 36000))
                              for x in (inclination, node, periapsis, anomaly)]
    p = a * (1 - e * e)
    r = p / (1 + e * math.cos(nu))
    mu = G * (star.mass + mass)
    v = math.sqrt(mu / p)
    planar_p = [r * math.cos(nu), r * math.sin(nu), 0]
    planar_v = [-v * math.sin(nu), v * (e + math.cos(nu)), 0]

    def rotate(q):
        x = math.cos(arg) * q[0] - math.sin(arg) * q[1]
        y = math.sin(arg) * q[0] + math.cos(arg) * q[1]
        y, z = y * math.cos(inc), y * math.sin(inc)
        return [math.cos(ascending) * x - math.sin(ascending) * y,
                math.sin(ascending) * x + math.cos(ascending) * y, z]

    return add(star.pos, rotate(planar_p)), add(star.vel, rotate(planar_v))


class CollisionError(ValueError):
    pass


class System:
    def __init__(self, bodies=None, epsilon=1e-7, stop_collision=True):
        self.bodies = list(bodies or [])
        self.epsilon = epsilon
        self.stop_collision = stop_collision
        self.days = 0.0

    def accelerations(self, positions):
        acc = [[0.0, 0.0, 0.0] for _ in self.bodies]
        for i in range(len(self.bodies)):
            for j in range(i + 1, len(self.bodies)):
                delta = sub(positions[j], positions[i])
                r2 = dot(delta, delta) + self.epsilon**2
                if r2 == 0:
                    raise ValueError("Corpos coincidentes: aumente a suavização ou separe-os.")
                factor = G / (r2 * math.sqrt(r2))
                for k in range(3):
                    value = delta[k] * factor
                    acc[i][k] += self.bodies[j].mass * value
                    acc[j][k] -= self.bodies[i].mass * value
        return acc

    def safe_step(self, requested):
        """Limita o passo a 3% das escalas locais de queda livre e cruzamento."""
        h = requested
        for i, a in enumerate(self.bodies):
            for b in self.bodies[i + 1:]:
                r = math.sqrt(dot(sub(b.pos, a.pos), sub(b.pos, a.pos)) + self.epsilon**2)
                if r == 0:
                    return 0.0
                h = min(h, 0.03 * math.sqrt(r**3 / (G * (a.mass + b.mass))))
                speed = norm(sub(b.vel, a.vel))
                if speed:
                    h = min(h, 0.03 * r / speed)
        return h

    def check_collisions(self, before, after):
        """Contato físico e interseção no segmento percorrido durante um passo."""
        if not self.stop_collision:
            return
        for i, a in enumerate(self.bodies):
            for j in range(i + 1, len(self.bodies)):
                b = self.bodies[j]
                radius = (a.radius_km + b.radius_km) / AU_KM
                d0 = sub(before[j], before[i])
                dd = sub(sub(after[j], after[i]), d0)
                d2 = dot(dd, dd)
                fraction = max(0.0, min(1.0, -dot(d0, dd) / d2)) if d2 else 0.0
                if norm(add(d0, mul(dd, fraction))) <= radius:
                    raise CollisionError("Contato físico entre %s e %s.\n"
                                         "A simulação parou antes deste passo.\n"
                                         "Edite as condições ou desmarque 'Pausar em colisão'."
                                         % (a.name, b.name))

    def step(self, h):
        if not math.isfinite(h) or h == 0:
            raise ValueError("Passo de tempo inválido.")
        old = [b.pos[:] for b in self.bodies]
        self.check_collisions(old, old)
        a0 = self.accelerations(old)
        new = [add(add(b.pos, mul(b.vel, h)), mul(a0[i], 0.5 * h * h))
               for i, b in enumerate(self.bodies)]
        self.check_collisions(old, new)
        a1 = self.accelerations(new)
        velocities = [add(b.vel, mul(add(a0[i], a1[i]), 0.5 * h))
                      for i, b in enumerate(self.bodies)]
        if not all(math.isfinite(x) for vec in new + velocities for x in vec):
            raise ValueError("Instabilidade numérica. Reduza o passo de integração.")
        for b, p, v in zip(self.bodies, new, velocities):
            b.pos, b.vel = p, v
        self.days += h

    def barycenter(self):
        total = sum(b.mass for b in self.bodies)
        if not total:
            return [0.0, 0.0, 0.0]
        return [sum(b.mass / total * b.pos[k] for b in self.bodies) for k in range(3)]

    def center_momentum(self):
        if not self.bodies:
            return
        total = sum(b.mass for b in self.bodies)
        center = self.barycenter()
        velocity = [sum(b.mass / total * b.vel[k] for b in self.bodies) for k in range(3)]
        for b in self.bodies:
            b.pos = sub(b.pos, center)
            b.vel = sub(b.vel, velocity)

    def energy(self):
        """Joules; potencial consistente com a suavização de Plummer."""
        value = sum(0.5 * b.mass * dot(b.vel, b.vel) for b in self.bodies)
        for i, a in enumerate(self.bodies):
            for b in self.bodies[i + 1:]:
                delta = sub(b.pos, a.pos)
                r = math.sqrt(dot(delta, delta) + self.epsilon**2)
                if r == 0:
                    raise ValueError("Energia indefinida para corpos coincidentes.")
                value -= G * a.mass * b.mass / r
        return value * (AU_KM * 1000.0 / DAY_S)**2


def scenario(name):
    sun = Body("Sol", SUN_KG, 696340, [0, 0, 0], [0, 0, 0],
               "Estrela", color="#ffcb62", visual_px=15)
    bodies = [sun]
    if name == "Sistema solar interno":
        rows = [("Mercúrio", 3.3011e23, 2439.7, .3871, .2056, 7.0, 48.3, 29.1, 45, "#b9b5ae", 4),
                ("Vênus", 4.8675e24, 6051.8, .7233, .0068, 3.39, 76.7, 54.9, 170, "#efbe81", 7),
                ("Terra", EARTH_KG, 6371, 1.0, .0167, 0, 0, 102.9, 250, "#57bfff", 7),
                ("Marte", 6.4171e23, 3389.5, 1.5237, .0934, 1.85, 49.6, 286.5, 25, "#f58265", 6)]
    elif name == "Órbitas em 3D":
        rows = [("Azul · 0°", EARTH_KG, 6371, 1, .02, 0, 0, 0, 0, "#59bfff", 8),
                ("Verde · 45°", EARTH_KG, 6371, 1.6, .25, 45, 60, 30, 120, "#74e5b2", 8),
                ("Rosa · 80°", EARTH_KG, 6371, 2.3, .4, 80, 140, 0, 230, "#ef91cf", 8)]
    elif name == "Estrelas binárias":
        sun.name = "Estrela A"
        other = Body("Estrela B", .7 * SUN_KG, 500000, [1.8, 0, 0], [0, 0, 0],
                     "Estrela", color="#a2c9ff", visual_px=12)
        other.vel[1] = math.sqrt(G * (sun.mass + other.mass) / 1.8)
        bodies.append(other)
        rows = [("Planeta externo", EARTH_KG, 6371, 7, .03, 18, 0, 0, 100, "#78ddc4", 7)]
        # Órbita circumbinária aproximada, centrada no baricentro das estrelas.
        system = System(bodies)
        system.center_momentum()
        center = Body("Centro", sun.mass + other.mass, 0, [0, 0, 0], [0, 0, 0])
        p, v = orbital_state(center, EARTH_KG, 7, .03, 18, 0, 0, 100)
        bodies.append(Body("Planeta externo", EARTH_KG, 6371, p, v, "Livre",
                           color="#78ddc4", visual_px=7))
        rows = []
    elif name == "Corpo livre em passagem":
        rows = [("Terra", EARTH_KG, 6371, 1, 0, 0, 0, 0, 150, "#57bfff", 7)]
        bodies.append(Body("Visitante", 1e18, 50, [-3, .7, .4],
                           mul([42, 0, -1], KM_S_TO_AU_DAY), color="#bf9aff", visual_px=5))
    elif name == "Vazio":
        return System([])
    else:
        raise ValueError("Exemplo desconhecido.")
    for title, mass, radius, a, e, inc, node, arg, nu, color, pixels in rows:
        p, v = orbital_state(sun, mass, a, e, inc, node, arg, nu)
        bodies.append(Body(title, mass, radius, p, v, "Em órbita", sun.uid, color, pixels))
    result = System(bodies)
    result.center_momentum()
    return result


def system_to_dict(system):
    return {"versao": 1, "tempo_dias": system.days,
            "suavizacao_ua": system.epsilon, "pausar_colisao": system.stop_collision,
            "corpos": [b.to_dict() for b in system.bodies]}


def system_from_dict(data):
    if not isinstance(data, dict) or data.get("versao") != 1:
        raise ValueError("Formato de cenário não reconhecido (versão 1 esperada).")
    records = data.get("corpos")
    if not isinstance(records, list) or len(records) > MAX_BODIES:
        raise ValueError("Lista de corpos inválida ou acima do limite.")
    bodies = [Body.from_dict(r) for r in records]
    validate_links(bodies)
    collision = data.get("pausar_colisao", True)
    if not isinstance(collision, bool):
        raise ValueError("pausar_colisao deve ser true ou false.")
    system = System(bodies, number(data.get("suavizacao_ua", 1e-7), "Suavização", 0, 1e6), collision)
    system.days = number(data.get("tempo_dias", 0), "Tempo", -1e12, 1e12)
    return system


def display_position(body, by_id):
    """Transformação somente visual. O pai de um corpo orbital é sempre estrela."""
    if body.kind == "Em órbita" and body.parent in by_id:
        parent = by_id[body.parent]
        anchor = mul(parent.pos, parent.distance_factor)
        return add(anchor, mul(sub(body.pos, parent.pos), body.distance_factor))
    return mul(body.pos, body.distance_factor)


@dataclass
class Camera:
    yaw: float = -25.0
    pitch: float = 58.0
    zoom: float = 170.0
    perspective: bool = True
    pan_x: float = 0.0
    pan_y: float = 0.0
    focal: float = 1100.0

    def project(self, point, center, width, height):
        x, y, z = sub(point, center)
        yaw, pitch = math.radians(self.yaw), math.radians(self.pitch)
        horizontal = math.cos(yaw) * x - math.sin(yaw) * y
        forward = math.sin(yaw) * x + math.cos(yaw) * y
        vertical = math.cos(pitch) * forward - math.sin(pitch) * z
        depth = math.sin(pitch) * forward + math.cos(pitch) * z
        scale = 1.0
        if self.perspective:
            denominator = self.focal + depth * self.zoom
            if denominator < self.focal * .06:
                return None
            scale = self.focal / denominator
        return (width / 2 + self.pan_x + horizontal * self.zoom * scale,
                height / 2 + self.pan_y - vertical * self.zoom * scale,
                depth, scale)


HELP = """SIMULADOR CELESTE 3D

INÍCIO RÁPIDO
Escolha um exemplo e clique em Carregar. Use Iniciar/Pausar; Passo avança um
passo físico. Reiniciar volta ao estado de referência. Novo corpo abre o cadastro.
Adicionar, editar, excluir ou mudar o referencial cria uma nova referência no
instante atual. Salvar grava o estado atual; Abrir o recupera pausado.

MOVIMENTO E TIPOS
Livre: posição e velocidade cartesianas definidas por você.
Estrela: mesmas leis físicas; pode servir de referência ao criar órbitas.
Em órbita: corpo associado a uma estrela; use coordenadas manuais ou o gerador
de órbita. A associação não prende o corpo: todos atraem todos, estrelas também
se movem. Um corpo associado pode escapar pela dinâmica do sistema.

GERADOR DE ÓRBITA
Escolha a estrela e informe semieixo maior a (UA), excentricidade e (0=circular),
inclinação, longitude do nó ascendente, argumento do periastro e anomalia
verdadeira. Os ângulos são em graus. Gerar preenche posição e velocidade
globais, somando o estado da estrela. A órbita é uma condição inicial de dois
corpos; outros corpos a perturbam. A estrela não é reajustada automaticamente
ao gerar. 'Referencial do baricentro' retira a translação do conjunto inteiro.

UNIDADES FÍSICAS
Massa: kg. Raio: km. Posição: UA. Velocidade: km/s. Tempo: dias.
1 UA = 149.597.870,7 km. Digite 5.9722e24 para a massa da Terra.
Pode usar vírgula decimal; não use separadores de milhares.

CÂMERA E ESCALAS
Arraste com o botão esquerdo para girar. Roda: zoom das distâncias.
Arraste com botão direito ou Shift+esquerdo para deslocar a vista.
Clique num corpo para selecioná-lo. Duplo clique na lista abre sua edição.
As vistas XY, XZ e YZ ajudam a inspecionar os três planos.
Escala de distância = pixels por UA no plano da câmera.
Tamanho global multiplica os raios visuais sem alterar as distâncias.
No modo Manual, cada corpo tem seu raio visual em pixels. No modo Físico,
raios em km usam a mesma escala das posições, vezes o multiplicador de tamanho.
Corpos abaixo de 1,5 px recebem um marcador mínimo para serem selecionáveis.
No cadastro, Fator distância visual expande/comprime a posição do corpo:
para 'Em órbita', em torno da estrela; para os demais, em torno da origem.
Essa distorção é apenas ilustrativa. Distâncias no painel são sempre físicas.
Trilhas guardam posições globais passadas. Alterar fatores visuais limpa trilhas.

PRECISÃO E LIMITES
Velocity Verlet, forças newtonianas em 3D e suavização de Plummer. O passo
máximo é reduzido automaticamente nas aproximações. Menor passo custa mais
processamento. Um passo variável perde a propriedade simplética estrita.
A velocidade em dias/s é uma meta: sob carga, o simulador desacelera.
Desvio de energia ajuda a avaliar erros; não certifica precisão sozinho.
Suavização evita singularidades e modifica a força em distâncias muito pequenas.
Colisões usam raios FÍSICOS. Por padrão o programa pausa antes de um passo
com contato; não há fusão ou ricochete. Sem pausa, corpos podem se atravessar.
O teste de colisão usa segmentos por passo, aproximação melhor com passos menores.
Sem relatividade, marés, rotação axial, atmosfera ou evolução estelar.
Exemplos são didáticos e não representam efemérides em uma data real.
Forças têm custo O(N²); 5 a 30 corpos são uma faixa prática, dependendo do PC.
Limite de cadastro: 150. Renderização 3D por projeção em Canvas, sem OpenGL;
discos sombreados representam esferas, com ordenação aproximada de profundidade.

EXECUÇÃO
Python 3.9+ com Tkinter/Tcl/Tk. Não há pip nem bibliotecas de terceiros.
No Windows, Tcl/Tk costuma acompanhar a instalação oficial do Python.
Algumas distribuições mínimas de Linux/macOS não incluem Tk: a interface exige
um Python já com Tkinter e uma sessão gráfica. Testes: --testar (sem interface).
"""


def launch_gui(initial=None):
    global tk, ttk, messagebox, filedialog, colorchooser
    try:
        import tkinter as tk
        from tkinter import ttk, messagebox, filedialog, colorchooser
    except ImportError as exc:
        print("Este Python não inclui Tkinter. Use um Python com Tcl/Tk incluído.\n%s" % exc, file=sys.stderr)
        return 1
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        print("Não foi possível abrir a interface gráfica. É necessário Tk e uma sessão com tela.\n%s" % exc,
              file=sys.stderr)
        return 1
    Application(root, initial)
    root.mainloop()
    return 0


class BodyDialog:
    def __init__(self, app, body=None):
        self.app, self.original = app, body
        self.result = None
        self.window = tk.Toplevel(app.root)
        self.window.title("Editar corpo" if body else "Novo corpo celeste")
        self.window.geometry("660x710")
        self.window.minsize(610, 650)
        self.window.transient(app.root)
        self.window.grab_set()
        self.vars = {}
        self.stars = [b for b in app.system.bodies if b.kind == "Estrela" and b is not body]
        self.star_labels = {"%s [%s]" % (b.name, b.uid[:6]): b.uid for b in self.stars}
        base = body or Body("Novo corpo", EARTH_KG, 6371, [1, 0, 0],
                            mul([0, 29.78, 0], KM_S_TO_AU_DAY))
        frame = ttk.Frame(self.window, padding=14)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)

        def entry(row, label, key, value):
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", pady=4)
            var = tk.StringVar(value=str(value))
            self.vars[key] = var
            ttk.Entry(frame, textvariable=var).grid(row=row, column=1, columnspan=2, sticky="ew", pady=4)

        entry(0, "Nome", "name", base.name)
        self.vars["kind"] = tk.StringVar(value=base.kind)
        ttk.Label(frame, text="Tipo").grid(row=1, column=0, sticky="w")
        types = ttk.Combobox(frame, textvariable=self.vars["kind"], values=KINDS, state="readonly")
        types.grid(row=1, column=1, columnspan=2, sticky="ew", pady=4)
        self.vars["parent"] = tk.StringVar(value=next((label for label, uid in self.star_labels.items()
                                                     if uid == base.parent), ""))
        ttk.Label(frame, text="Estrela associada").grid(row=2, column=0, sticky="w")
        self.star_box = ttk.Combobox(frame, textvariable=self.vars["parent"],
                                    values=list(self.star_labels), state="readonly")
        self.star_box.grid(row=2, column=1, columnspan=2, sticky="ew", pady=4)
        entry(3, "Massa (kg)", "mass", "%g" % base.mass)
        entry(4, "Raio físico (km)", "radius", "%g" % base.radius_km)
        entry(5, "Raio visual manual (px)", "pixels", "%g" % base.visual_px)
        entry(6, "Fator distância visual", "distance", "%g" % base.distance_factor)
        self.vars["color"] = tk.StringVar(value=base.color)
        ttk.Label(frame, text="Cor (#RRGGBB)").grid(row=7, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.vars["color"]).grid(row=7, column=1, sticky="ew")
        ttk.Button(frame, text="Escolher", command=self.choose_color).grid(row=7, column=2, padx=(6, 0))
        tabs = ttk.Notebook(frame)
        tabs.grid(row=8, column=0, columnspan=3, sticky="nsew", pady=12)
        frame.rowconfigure(8, weight=1)
        cart = ttk.Frame(tabs, padding=12)
        orbit = ttk.Frame(tabs, padding=12)
        tabs.add(cart, text="Posição e velocidade")
        tabs.add(orbit, text="Gerador de órbita")
        ttk.Label(cart, text="Coordenadas globais. Z permite movimento fora do plano.").grid(
            row=0, column=0, columnspan=4, sticky="w", pady=(0, 12))
        for k, axis in enumerate("XYZ"):
            ttk.Label(cart, text=axis).grid(row=1, column=k + 1)
        for row, (label, prefix, vec) in enumerate((
                ("Posição (UA)", "p", base.pos),
                ("Velocidade (km/s)", "v", mul(base.vel, 1 / KM_S_TO_AU_DAY))), 2):
            ttk.Label(cart, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10))
            for k in range(3):
                var = tk.StringVar(value="%.12g" % vec[k])
                self.vars[prefix + str(k)] = var
                ttk.Entry(cart, textvariable=var, width=15).grid(row=row, column=k + 1, padx=3, pady=5)
        ttk.Label(cart, text="Exemplo: Terra ≈ 5.9722e24 kg; 1 UA; 29,78 km/s.\n"
                            "Livre e Estrela também se movem e sofrem gravidade.",
                  wraplength=535).grid(row=4, column=0, columnspan=4, sticky="w", pady=18)
        orbital = [("a", "Semieixo maior (UA)", 1), ("e", "Excentricidade [0; 0,99]", 0),
                   ("i", "Inclinação (°)", 0), ("node", "Longitude do nó (°)", 0),
                   ("arg", "Argumento do periastro (°)", 0), ("nu", "Anomalia verdadeira (°)", 0)]
        for row, (key, label, default) in enumerate(orbital):
            ttk.Label(orbit, text=label).grid(row=row, column=0, sticky="w", pady=3)
            self.vars[key] = tk.StringVar(value=str(default))
            ttk.Entry(orbit, textvariable=self.vars[key], width=20).grid(row=row, column=1, padx=14)
        ttk.Button(orbit, text="Gerar posição e velocidade", command=lambda: self.generate(tabs, cart)).grid(
            row=6, column=0, columnspan=2, sticky="ew", pady=10)
        ttk.Label(frame, text="O gerador altera os campos; clique em Salvar corpo para aplicar.\n"
                              "Escalas visuais não alteram a força gravitacional.",
                  wraplength=600).grid(row=9, column=0, columnspan=3, sticky="w", pady=5)
        buttons = ttk.Frame(frame)
        buttons.grid(row=10, column=0, columnspan=3, sticky="e", pady=(10, 0))
        ttk.Button(buttons, text="Cancelar", command=self.window.destroy).pack(side="left", padx=6)
        ttk.Button(buttons, text="Salvar corpo", command=self.save).pack(side="left")
        self.window.bind("<Escape>", lambda event: self.window.destroy())
        types.bind("<<ComboboxSelected>>", lambda event: self.update_kind())
        self.update_kind()

    def update_kind(self):
        orbital = self.vars["kind"].get() == "Em órbita"
        self.star_box.configure(state="readonly" if orbital else "disabled")
        if orbital and not self.vars["parent"].get() and self.star_labels:
            self.vars["parent"].set(next(iter(self.star_labels)))

    def choose_color(self):
        current = self.vars["color"].get()
        if not re.fullmatch(r"#[0-9a-fA-F]{6}", current):
            current = "#8dccff"
        result = colorchooser.askcolor(current, parent=self.window)[1]
        if result:
            self.vars["color"].set(result)

    def generate(self, tabs, cart):
        try:
            if not self.stars:
                raise ValueError("Cadastre uma estrela antes de gerar uma órbita.")
            self.vars["kind"].set("Em órbita")
            self.update_kind()
            uid = self.star_labels.get(self.vars["parent"].get())
            star = next((b for b in self.stars if b.uid == uid), None)
            if star is None:
                raise ValueError("Escolha a estrela associada.")
            mass = number(self.vars["mass"].get(), "Massa", 1e-6, 1e36)
            p, v = orbital_state(star, mass, *[self.vars[k].get() for k in ("a", "e", "i", "node", "arg", "nu")])
            radius = number(self.vars["radius"].get(), "Raio", 0, 1e12)
            peri = number(self.vars["a"].get(), "a") * (1 - number(self.vars["e"].get(), "e"))
            if peri <= (radius + star.radius_km) / AU_KM:
                raise ValueError("O periastro dessa órbita atravessa o raio físico da estrela.")
            for prefix, values in (("p", p), ("v", mul(v, 1 / KM_S_TO_AU_DAY))):
                for k in range(3):
                    self.vars[prefix + str(k)].set("%.12g" % values[k])
            tabs.select(cart)
        except ValueError as exc:
            messagebox.showerror("Verifique os dados", str(exc), parent=self.window)

    def save(self):
        try:
            get = lambda key: self.vars[key].get()
            kind = get("kind")
            body = Body.from_dict({
                "id": self.original.uid if self.original else uuid.uuid4().hex,
                "nome": get("name"), "tipo": kind,
                "estrela_id": self.star_labels.get(get("parent"), "") if kind == "Em órbita" else "",
                "massa_kg": get("mass"), "raio_km": get("radius"),
                "posicao_ua": [get("p" + str(i)) for i in range(3)],
                "velocidade_km_s": [get("v" + str(i)) for i in range(3)],
                "cor": get("color"), "raio_visual_px": get("pixels"),
                "fator_distancia_visual": get("distance")})
            bodies = [body if b.uid == body.uid else b for b in self.app.system.bodies]
            if self.original is None:
                bodies.append(body)
            validate_links(bodies)
            self.result = body
            self.window.destroy()
        except ValueError as exc:
            messagebox.showerror("Verifique os dados", str(exc), parent=self.window)


class Application:
    def __init__(self, root, initial=None):
        self.root = root
        root.title("Simulador Celeste 3D · gravitação de múltiplos corpos")
        root.geometry("1380x900")
        root.minsize(1030, 690)
        self.system = initial or scenario("Sistema solar interno")
        self.camera = Camera()
        self.running = False
        self.selected = ""
        self.trails = {}
        self.hits = []
        self.drag = None
        self.actual_speed = 0.0
        self.last_ui = 0.0
        self.last_trail_day = -math.inf
        self.last_clock = time.perf_counter()
        self.vars = {}
        self.dirty = False
        style = ttk.Style(root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("TButton", padding=(9, 5))
        style.configure("TLabel", font=("Segoe UI", 10))
        style.configure("Title.TLabel", font=("Segoe UI", 16, "bold"))
        style.configure("Treeview", rowheight=25)
        self.build_ui()
        self.set_reference()
        self.refresh_list()
        self.fit()
        root.after_idle(self.fit)
        root.protocol("WM_DELETE_WINDOW", self.close)
        root.after(30, self.tick)

    def var(self, key, value, kind=None):
        self.vars[key] = (kind or tk.StringVar)(value=value)
        return self.vars[key]

    def build_ui(self):
        top = ttk.Frame(self.root, padding=(12, 10))
        top.pack(fill="x")
        ttk.Label(top, text="SIMULADOR CELESTE", style="Title.TLabel").pack(side="left")
        for label, command in (("Ajuda", self.help), ("Salvar", self.save), ("Abrir", self.open)):
            ttk.Button(top, text=label, command=command).pack(side="right", padx=3)
        bar = ttk.Frame(self.root, padding=(12, 0, 12, 9))
        bar.pack(fill="x")
        self.run_button = ttk.Button(bar, text="▶ Iniciar", command=self.toggle)
        self.run_button.pack(side="left")
        ttk.Button(bar, text="Passo", command=self.single_step).pack(side="left", padx=4)
        ttk.Button(bar, text="Reiniciar", command=self.reset).pack(side="left")
        ttk.Label(bar, text="   dias/s").pack(side="left")
        ttk.Entry(bar, textvariable=self.var("speed", "10"), width=8).pack(side="left", padx=5)
        ttk.Label(bar, text="Passo máx. (dias)").pack(side="left", padx=(8, 0))
        ttk.Entry(bar, textvariable=self.var("step", "0.05"), width=8).pack(side="left", padx=5)
        self.time_label = ttk.Label(bar, text="")
        self.time_label.pack(side="right")
        panes = ttk.Panedwindow(self.root, orient="horizontal")
        panes.pack(fill="both", expand=True, padx=12)
        left = ttk.Frame(panes, width=335)
        right = ttk.Frame(panes)
        panes.add(left, weight=0)
        panes.add(right, weight=1)
        # Sidebar rolável: acessível também em telas menores e com escala de DPI.
        sidebar = tk.Canvas(left, width=330, highlightthickness=0, background="#e8e8e8")
        scrollbar = ttk.Scrollbar(left, orient="vertical", command=sidebar.yview)
        sidebar.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        sidebar.pack(fill="both", expand=True)
        controls = ttk.Frame(sidebar, padding=(4, 4, 10, 12))
        window_id = sidebar.create_window((0, 0), window=controls, anchor="nw")
        controls.bind("<Configure>", lambda e: sidebar.configure(scrollregion=sidebar.bbox("all")))
        sidebar.bind("<Configure>", lambda e: sidebar.itemconfigure(window_id, width=e.width))
        examples = ttk.LabelFrame(controls, text="Cenários", padding=8)
        examples.pack(fill="x", pady=(0, 8))
        self.example = ttk.Combobox(examples, state="readonly", values=(
            "Sistema solar interno", "Órbitas em 3D", "Estrelas binárias", "Corpo livre em passagem", "Vazio"))
        self.example.set("Sistema solar interno")
        self.example.pack(fill="x")
        ttk.Button(examples, text="Carregar exemplo", command=self.load_example).pack(fill="x", pady=(6, 0))
        bodies = ttk.LabelFrame(controls, text="Corpos", padding=8)
        bodies.pack(fill="x", pady=(0, 8))
        self.tree = ttk.Treeview(bodies, columns=("kind",), show="tree headings", height=6, selectmode="browse")
        self.tree.heading("#0", text="Nome")
        self.tree.heading("kind", text="Tipo")
        self.tree.column("#0", width=145, minwidth=80)
        self.tree.column("kind", width=85, minwidth=65)
        tree_scroll = ttk.Scrollbar(bodies, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=tree_scroll.set)
        tree_scroll.pack(side="right", fill="y")
        self.tree.pack(fill="x")
        self.tree.bind("<<TreeviewSelect>>", self.on_select)
        self.tree.bind("<Double-1>", lambda e: self.edit_body())
        actions = ttk.Frame(bodies)
        actions.pack(fill="x", pady=(7, 0))
        for label, callback in (("+ Novo", lambda: self.edit_body(True)), ("Editar", self.edit_body), ("Excluir", self.delete_body)):
            ttk.Button(actions, text=label, command=callback).pack(side="left", expand=True, fill="x")
        visual = ttk.LabelFrame(controls, text="Visualização · não altera a física", padding=8)
        visual.pack(fill="x", pady=(0, 8))
        ttk.Label(visual, text="Distância: pixels / UA").pack(anchor="w")
        row = ttk.Frame(visual)
        row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var("zoom", "170"), width=12).pack(side="left")
        ttk.Button(row, text="Aplicar", command=self.apply_visual).pack(side="left", padx=5)
        ttk.Button(row, text="Enquadrar", command=self.fit).pack(side="left")
        ttk.Label(visual, text="Tamanho global: multiplicador").pack(anchor="w", pady=(8, 0))
        row = ttk.Frame(visual)
        row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var("size", "1"), width=12).pack(side="left")
        self.size_mode = ttk.Combobox(row, textvariable=self.var("size_mode", "Manual"),
                                     values=("Manual", "Físico"), state="readonly", width=12)
        self.size_mode.pack(side="left", padx=5)
        self.size_mode.bind("<<ComboboxSelected>>", lambda e: self.apply_visual())
        ttk.Label(visual, text="Use Aplicar após editar escala ou tamanho.", wraplength=285).pack(anchor="w", pady=4)
        self.size_multiplier = 1.0
        self.physical_size = False
        self.focus_ids = {"Origem": "", "Baricentro físico": "@bary"}
        ttk.Label(visual, text="Centro da câmera").pack(anchor="w", pady=(6, 0))
        self.focus_box = ttk.Combobox(visual, textvariable=self.var("focus", "Origem"), state="readonly")
        self.focus_box.pack(fill="x")
        self.focus_box.bind("<<ComboboxSelected>>", lambda e: self.focus_changed())
        ttk.Checkbutton(visual, text="Perspectiva", variable=self.var("perspective", True, tk.BooleanVar),
                        command=self.set_projection).pack(anchor="w", pady=(5, 0))
        for key, label, value in (("trails", "Trilhas", True), ("names", "Nomes", True), ("grid", "Grade e eixos", True)):
            ttk.Checkbutton(visual, text=label, variable=self.var(key, value, tk.BooleanVar)).pack(anchor="w")
        views = ttk.Frame(visual)
        views.pack(fill="x", pady=5)
        for name, yaw, pitch in (("3D", -25, 58), ("XY", 0, 0), ("XZ", 0, 90), ("YZ", -90, 90)):
            ttk.Button(views, text=name, width=5, command=lambda y=yaw, p=pitch: self.set_view(y, p)).pack(side="left")
        ttk.Button(visual, text="Limpar trilhas", command=self.clear_trails).pack(fill="x")
        physical = ttk.LabelFrame(controls, text="Integração e referencial", padding=8)
        physical.pack(fill="x", pady=(0, 8))
        ttk.Label(physical, text="Suavização gravitacional (UA)").pack(anchor="w")
        row = ttk.Frame(physical)
        row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var("epsilon", "%g" % self.system.epsilon), width=14).pack(side="left")
        ttk.Button(row, text="Aplicar", command=self.apply_physics).pack(side="left", padx=5)
        ttk.Checkbutton(physical, text="Pausar em colisão física", variable=self.var(
            "collision", self.system.stop_collision, tk.BooleanVar), command=self.apply_physics).pack(anchor="w", pady=7)
        ttk.Button(physical, text="Referencial do baricentro", command=self.recenter).pack(fill="x")
        self.info = ttk.Label(controls, text="", justify="left", wraplength=290)
        self.info.pack(fill="x", pady=8)
        self.canvas = tk.Canvas(right, background="#080f20", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<ButtonPress-1>", self.mouse_down)
        self.canvas.bind("<B1-Motion>", self.mouse_move)
        self.canvas.bind("<ButtonRelease-1>", self.mouse_up)
        self.canvas.bind("<ButtonPress-3>", self.mouse_down)
        self.canvas.bind("<B3-Motion>", self.mouse_move)
        self.canvas.bind("<ButtonRelease-3>", self.mouse_up)
        self.canvas.bind("<MouseWheel>", self.mouse_wheel)
        self.canvas.bind("<Button-4>", lambda e: self.zoom_by(1.12))
        self.canvas.bind("<Button-5>", lambda e: self.zoom_by(1 / 1.12))
        self.status = ttk.Label(self.root, text="", anchor="w", padding=(12, 7))
        self.status.pack(fill="x")

    def body(self):
        return next((b for b in self.system.bodies if b.uid == self.selected), None)

    def pause(self):
        self.running = False
        self.run_button.configure(text="▶ Iniciar")

    def controls(self):
        return (number(self.vars["speed"].get(), "Velocidade (dias/s)", 1e-8, 1e8),
                number(self.vars["step"].get(), "Passo máximo", 1e-10, 1e6))

    def toggle(self):
        if self.running:
            self.pause()
            return
        try:
            self.controls()
            if not self.system.bodies:
                raise ValueError("Adicione pelo menos um corpo.")
            self.running = True
            self.last_clock = time.perf_counter()
            self.run_button.configure(text="Ⅱ Pausar")
        except ValueError as exc:
            messagebox.showerror("Não foi possível iniciar", str(exc), parent=self.root)

    def single_step(self):
        self.pause()
        try:
            _, step = self.controls()
            if not self.system.bodies:
                raise ValueError("Adicione pelo menos um corpo.")
            h = self.system.safe_step(step)
            if h < 1e-12:
                raise ValueError("Encontro muito próximo: separe os corpos ou aumente a suavização.")
            self.system.step(h)
            self.dirty = True
            self.sample_trails(force=True)
        except ValueError as exc:
            messagebox.showerror("Simulação pausada", str(exc), parent=self.root)

    def clear_trails(self):
        self.trails = {b.uid: deque(maxlen=850) for b in self.system.bodies}
        self.last_trail_day = -math.inf
        self.sample_trails(force=True)

    def set_reference(self):
        self.reference = system_to_dict(self.system)
        self.energy_reference = self.get_energy()
        self.clear_trails()

    def get_energy(self):
        try:
            return self.system.energy()
        except ValueError:
            return None

    def reset(self):
        self.pause()
        self.system = system_from_dict(self.reference)
        self.vars["epsilon"].set("%g" % self.system.epsilon)
        self.vars["collision"].set(self.system.stop_collision)
        self.set_reference()
        self.refresh_list()
        self.dirty = True

    def confirm_discard(self):
        self.pause()
        return not self.dirty or messagebox.askyesno(
            "Substituir cenário?", "Há alterações ou evolução não salvas.\nContinuar sem salvar?", parent=self.root)

    def load_example(self):
        if not self.confirm_discard():
            return
        self.system = scenario(self.example.get())
        self.selected = ""
        self.vars["focus"].set("Origem")
        self.vars["epsilon"].set("%g" % self.system.epsilon)
        self.vars["collision"].set(self.system.stop_collision)
        self.set_reference()
        self.refresh_list()
        self.fit()
        self.dirty = False

    def refresh_list(self):
        self.tree.delete(*self.tree.get_children())
        self.focus_ids = {"Origem": "", "Baricentro físico": "@bary"}
        for b in self.system.bodies:
            self.tree.insert("", "end", iid=b.uid, text=b.name, values=(b.kind,))
            self.focus_ids["%s [%s]" % (b.name, b.uid[:6])] = b.uid
        self.focus_box.configure(values=list(self.focus_ids))
        if self.vars["focus"].get() not in self.focus_ids:
            self.vars["focus"].set("Origem")
        if self.body():
            self.tree.selection_set(self.selected)
        elif self.system.bodies:
            self.selected = self.system.bodies[0].uid
            self.tree.selection_set(self.selected)
        else:
            self.selected = ""

    def on_select(self, event=None):
        selection = self.tree.selection()
        if selection:
            self.selected = selection[0]

    def edit_body(self, new=False):
        self.pause()
        old = None if new else self.body()
        if not new and old is None:
            return
        if new and len(self.system.bodies) >= MAX_BODIES:
            messagebox.showerror("Limite", "Máximo de %d corpos." % MAX_BODIES, parent=self.root)
            return
        dialog = BodyDialog(self, old)
        self.root.wait_window(dialog.window)
        if dialog.result:
            body = dialog.result
            if old:
                self.system.bodies[self.system.bodies.index(old)] = body
            else:
                self.system.bodies.append(body)
            self.selected = body.uid
            self.set_reference()
            self.refresh_list()
            self.dirty = True

    def delete_body(self):
        b = self.body()
        if b is None:
            return
        self.pause()
        dependents = [c for c in self.system.bodies if c.parent == b.uid]
        detail = "\nOs corpos associados passarão a ser livres." if dependents else ""
        if not messagebox.askyesno("Excluir", "Excluir %s?%s" % (b.name, detail), parent=self.root):
            return
        self.system.bodies.remove(b)
        for child in dependents:
            child.kind, child.parent = "Livre", ""
        self.set_reference()
        self.refresh_list()
        self.dirty = True

    def apply_visual(self):
        try:
            zoom = number(self.vars["zoom"].get(), "Pixels por UA", 1e-8, 1e12)
            size = number(self.vars["size"].get(), "Multiplicador do tamanho", 1e-4, 1e8)
            self.camera.zoom, self.size_multiplier = zoom, size
            self.physical_size = self.vars["size_mode"].get() == "Físico"
            self.dirty = True
        except ValueError as exc:
            messagebox.showerror("Escala inválida", str(exc), parent=self.root)

    def apply_physics(self):
        self.pause()
        try:
            epsilon = number(self.vars["epsilon"].get(), "Suavização", 0, 1e6)
            self.system.epsilon = epsilon
            self.system.stop_collision = self.vars["collision"].get()
            self.set_reference()
            self.dirty = True
        except ValueError as exc:
            self.vars["collision"].set(self.system.stop_collision)
            messagebox.showerror("Valor inválido", str(exc), parent=self.root)

    def recenter(self):
        self.pause()
        self.system.center_momentum()
        self.set_reference()
        self.dirty = True

    def center(self):
        uid = self.focus_ids.get(self.vars["focus"].get(), "")
        if uid == "@bary":
            return self.system.barycenter()
        by_id = {b.uid: b for b in self.system.bodies}
        return display_position(by_id[uid], by_id) if uid in by_id else [0, 0, 0]

    def focus_changed(self):
        self.camera.pan_x = self.camera.pan_y = 0
        self.dirty = True

    def set_projection(self):
        self.camera.perspective = self.vars["perspective"].get()
        self.dirty = True

    def set_view(self, yaw, pitch):
        self.camera.yaw, self.camera.pitch = yaw, pitch
        self.camera.pan_x = self.camera.pan_y = 0
        self.dirty = True

    def fit(self):
        center = self.center()
        by_id = {b.uid: b for b in self.system.bodies}
        extent = max([norm(sub(display_position(b, by_id), center)) for b in self.system.bodies] or [2])
        width, height = max(500, self.canvas.winfo_width()), max(400, self.canvas.winfo_height())
        self.camera.zoom = max(1e-8, min(1e12, min(width, height, 1100) * .28 / max(extent, 1e-9)))
        self.camera.pan_x = self.camera.pan_y = 0
        self.vars["zoom"].set("%.6g" % self.camera.zoom)

    def mouse_down(self, event):
        self.drag = [event.x, event.y, event.x, event.y, event.num == 3 or bool(event.state & 0x0001), False]

    def mouse_move(self, event):
        if self.drag is None:
            return
        x, y, start_x, start_y, pan, moved = self.drag
        dx, dy = event.x - x, event.y - y
        moved = moved or abs(event.x - start_x) + abs(event.y - start_y) > 4
        if moved:
            if pan:
                self.camera.pan_x += dx
                self.camera.pan_y += dy
            else:
                self.camera.yaw = (self.camera.yaw + dx * .4) % 360
                self.camera.pitch = max(-90, min(90, self.camera.pitch + dy * .4))
            self.dirty = True
        self.drag = [event.x, event.y, start_x, start_y, pan, moved]

    def mouse_up(self, event):
        if self.drag is not None and not self.drag[5] and not self.drag[4]:
            for x, y, r, uid in reversed(self.hits):
                if math.hypot(event.x - x, event.y - y) < max(8, r):
                    self.selected = uid
                    self.tree.selection_set(uid)
                    self.tree.see(uid)
                    break
        self.drag = None

    def zoom_by(self, factor):
        self.camera.zoom = max(1e-8, min(1e12, self.camera.zoom * factor))
        self.vars["zoom"].set("%.6g" % self.camera.zoom)
        self.dirty = True

    def mouse_wheel(self, event):
        if event.delta:
            self.zoom_by(1.12 if event.delta > 0 else 1 / 1.12)

    def sample_trails(self, force=False):
        if not force and abs(self.system.days - self.last_trail_day) < .01:
            return
        by_id = {b.uid: b for b in self.system.bodies}
        for b in self.system.bodies:
            self.trails.setdefault(b.uid, deque(maxlen=850)).append(display_position(b, by_id))
        self.last_trail_day = self.system.days

    def tick(self):
        clock = time.perf_counter()
        wall = min(.15, clock - self.last_clock)
        self.last_clock = clock
        advanced = 0.0
        if self.running:
            try:
                speed, max_step = self.controls()
                target = speed * wall
                deadline = clock + .018
                count = 0
                while advanced < target and time.perf_counter() < deadline and count < 400:
                    h = self.system.safe_step(min(max_step, target - advanced))
                    if h < 1e-12:
                        if target - advanced < 1e-12:
                            break
                        raise ValueError("Passo necessário muito pequeno. Separe os corpos ou ajuste a suavização.")
                    self.system.step(h)
                    advanced += h
                    count += 1
                if advanced:
                    self.sample_trails()
                    self.dirty = True
                self.actual_speed = .85 * self.actual_speed + .15 * advanced / max(wall, 1e-9)
            except (ValueError, OverflowError) as exc:
                self.pause()
                messagebox.showerror("Simulação pausada", str(exc), parent=self.root)
        self.render()
        if clock - self.last_ui > .25:
            self.update_info()
            self.last_ui = clock
        self.root.after(30, self.tick)

    @staticmethod
    def tint(color, factor):
        rgb = [int(color[i:i + 2], 16) for i in (1, 3, 5)]
        return "#" + "".join("%02x" % max(0, min(255, int(v * factor))) for v in rgb)

    def render(self):
        c = self.canvas
        c.delete("all")
        width, height = c.winfo_width(), c.winfo_height()
        if width < 10 or height < 10:
            return
        center = self.center()
        project = lambda p: self.camera.project(p, center, width, height)
        if self.vars["grid"].get():
            raw = 90 / self.camera.zoom
            base = 10**math.floor(math.log10(raw))
            spacing = base * next((f for f in (1, 2, 5, 10) if f * base >= raw), 10)
            gx = round(center[0] / spacing) * spacing
            gy = round(center[1] / spacing) * spacing
            extent = 10 * spacing
            for j in range(-10, 11):
                self.line3d([gx + j * spacing, gy - extent, 0], [gx + j * spacing, gy + extent, 0], project, "#17263e")
                self.line3d([gx - extent, gy + j * spacing, 0], [gx + extent, gy + j * spacing, 0], project, "#17263e")
            for axis, color in enumerate(("#b45d70", "#54a58d", "#658ed1")):
                end = [0, 0, 0]
                end[axis] = spacing * 3
                self.line3d([0, 0, 0], end, project, color, 2)
                p = project(end)
                if p:
                    c.create_text(p[0] + 8, p[1], text="XYZ"[axis], fill=color, font=("Segoe UI", 10, "bold"))
            c.create_text(18, height - 24, anchor="w", text="Grade XY: %g UA   ·   eixos físicos" % spacing,
                          fill="#758aa8", font=("Segoe UI", 10))
        by_id = {b.uid: b for b in self.system.bodies}
        if self.vars["trails"].get():
            for b in self.system.bodies:
                points = list(self.trails.get(b.uid, []))
                stride = max(1, len(points) // 300)
                segment = []
                for point in points[::stride]:
                    p = project(point)
                    if p is None or abs(p[0]) > 1e6 or abs(p[1]) > 1e6:
                        if len(segment) >= 4:
                            c.create_line(*segment, fill=self.tint(b.color, .5), width=1)
                        segment = []
                    else:
                        segment.extend(p[:2])
                if len(segment) >= 4:
                    c.create_line(*segment, fill=self.tint(b.color, .5), width=1)
        projected = []
        for b in self.system.bodies:
            p = project(display_position(b, by_id))
            if p:
                projected.append((p[2], b, p))
        projected.sort(key=lambda item: item[0], reverse=True)
        self.hits = []
        for _, b, p in projected:
            x, y, _, perspective_scale = p
            basic = b.radius_km / AU_KM * self.camera.zoom if self.physical_size else b.visual_px
            radius = max(1.5, basic * self.size_multiplier * perspective_scale)
            # Culling impede coordenadas extremas de chegarem ao Tk.
            if x + radius < 0 or x - radius > width or y + radius < 0 or y - radius > height:
                continue
            radius = min(radius, 4 * max(width, height))
            if abs(x) > 1e6 or abs(y) > 1e6:
                continue
            if b.kind == "Estrela":
                c.create_oval(x - radius * 1.3, y - radius * 1.3, x + radius * 1.3, y + radius * 1.3,
                              fill=self.tint(b.color, .18), outline="")
            c.create_oval(x - radius, y - radius, x + radius, y + radius, fill=self.tint(b.color, .63), outline="")
            c.create_oval(x - radius * .85, y - radius * .87, x + radius * .6, y + radius * .6, fill=b.color, outline="")
            if radius >= 4:
                c.create_oval(x - radius * .58, y - radius * .58, x - radius * .13, y - radius * .13,
                              fill=self.tint(b.color, 1.2), outline="")
            if b.uid == self.selected:
                c.create_oval(x - radius - 5, y - radius - 5, x + radius + 5, y + radius + 5,
                              outline="#f3f7ff", width=1, dash=(3, 3))
            if self.vars["names"].get():
                c.create_text(x + radius + 8, y - radius - 3, anchor="sw", text=b.name,
                              fill="#dde8fa", font=("Segoe UI", 10))
            self.hits.append((x, y, radius, b.uid))
        c.create_text(18, 20, anchor="nw", text="GRAVITAÇÃO EM 3D", fill="#d2e2fa", font=("Segoe UI", 12, "bold"))
        c.create_text(18, 46, anchor="nw", text="Arraste: girar   ·   Roda: zoom   ·   Botão direito: mover",
                      fill="#839bbd", font=("Segoe UI", 10))
        if any(b.distance_factor != 1 for b in self.system.bodies):
            c.create_text(18, 69, anchor="nw", text="DISTÂNCIAS VISUAIS PERSONALIZADAS", fill="#f7c56f", font=("Segoe UI", 10))
        if not self.system.bodies:
            c.create_text(width / 2, height / 2, text="Adicione uma estrela ou um corpo livre para começar.",
                          fill="#a9c2e8", font=("Segoe UI", 14))

    def line3d(self, a, b, project, color, width=1):
        pa, pb = project(a), project(b)
        if pa and pb and max(abs(v) for v in (*pa[:2], *pb[:2])) < 1e6:
            self.canvas.create_line(*pa[:2], *pb[:2], fill=color, width=width)

    def update_info(self):
        self.time_label.configure(text="t = %.4f dias   |   %.3f anos" % (self.system.days, self.system.days / 365.25))
        b = self.body()
        if b:
            lines = [b.name + " · " + b.kind, "Massa: %.5g kg  |  Raio: %g km" % (b.mass, b.radius_km),
                     "Posição (UA):\n  x=%.6g  y=%.6g  z=%.6g" % tuple(b.pos),
                     "Velocidade: %.6g km/s" % (norm(b.vel) / KM_S_TO_AU_DAY)]
            parent = next((x for x in self.system.bodies if x.uid == b.parent), None)
            if parent:
                distance = norm(sub(b.pos, parent.pos))
                lines.append("Até %s: %.6g UA\n  (%.6g km)" % (parent.name, distance, distance * AU_KM))
            others = [x for x in self.system.bodies if x.uid != b.uid]
            if others:
                nearest = min(others, key=lambda x: norm(sub(x.pos, b.pos)))
                lines.append("Mais próximo: %s\n  %.6g UA (centro a centro)" % (nearest.name, norm(sub(nearest.pos, b.pos))))
            self.info.configure(text="\n\n".join(lines))
        else:
            self.info.configure(text="Selecione um corpo para ver os dados físicos.")
        energy = self.get_energy()
        if energy is not None and self.energy_reference is not None and abs(self.energy_reference) > 1e-30:
            error = abs((energy - self.energy_reference) / self.energy_reference)
            energy_text = "Desvio de energia: %.3g%%" % (100 * error)
        else:
            energy_text = "Desvio de energia: indisponível"
        self.status.configure(text="%s   |   %d corpos   |   %.3g dias/s efetivos   |   %s" % (
            "Simulando" if self.running else "Pausado", len(self.system.bodies),
            self.actual_speed if self.running else 0, energy_text))

    def camera_dict(self):
        return {"yaw": self.camera.yaw, "pitch": self.camera.pitch, "zoom": self.camera.zoom,
                "pan_x": self.camera.pan_x, "pan_y": self.camera.pan_y,
                "perspectiva": self.camera.perspective, "tamanho": self.size_multiplier,
                "tamanho_fisico": self.physical_size,
                "foco": self.focus_ids.get(self.vars["focus"].get(), ""),
                "trilhas": self.vars["trails"].get(), "nomes": self.vars["names"].get(),
                "grade": self.vars["grid"].get()}

    def save(self):
        self.pause()
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".json",
                                          filetypes=[("Cenário JSON", "*.json")], initialfile="meu_universo.json")
        if not path:
            return
        try:
            speed, step = self.controls()
            data = system_to_dict(self.system)
            data["camera"] = self.camera_dict()
            data["execucao"] = {"dias_por_segundo": speed, "passo_max_dias": step}
            # Substituição atômica para evitar um arquivo final parcialmente escrito.
            import tempfile
            import os
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=str(Path(path).parent),
                                                 delete=False, suffix=".tmp") as stream:
                    temporary = stream.name
                    json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
                    stream.write("\n")
                os.replace(temporary, path)
            finally:
                if temporary and os.path.exists(temporary):
                    os.unlink(temporary)
            self.dirty = False
        except (OSError, ValueError) as exc:
            messagebox.showerror("Não foi possível salvar", str(exc), parent=self.root)

    def open(self):
        self.pause()
        path = filedialog.askopenfilename(parent=self.root, filetypes=[("Cenário JSON", "*.json")])
        if not path:
            return
        try:
            data = read_json(path)
            system = system_from_dict(data)
            camera, extras = parse_view(data)
            if not self.confirm_discard():
                return
            self.system, self.camera = system, camera
            self.selected = ""
            self.vars["epsilon"].set("%g" % system.epsilon)
            self.vars["collision"].set(system.stop_collision)
            self.size_multiplier, self.physical_size = extras["size"], extras["physical"]
            self.vars["size"].set("%g" % self.size_multiplier)
            self.vars["size_mode"].set("Físico" if self.physical_size else "Manual")
            self.vars["zoom"].set("%.6g" % camera.zoom)
            self.vars["perspective"].set(camera.perspective)
            self.vars["speed"].set("%g" % extras["speed"])
            self.vars["step"].set("%g" % extras["step"])
            for key in ("trails", "names", "grid"):
                self.vars[key].set(extras[key])
            self.set_reference()
            self.refresh_list()
            self.vars["focus"].set(next((k for k, uid in self.focus_ids.items() if uid == extras["focus"]), "Origem"))
            if "camera" not in data:
                self.fit()
            self.dirty = False
        except (OSError, ValueError, TypeError) as exc:
            messagebox.showerror("Não foi possível abrir", str(exc), parent=self.root)

    def help(self):
        self.pause()
        window = tk.Toplevel(self.root)
        window.title("Guia do simulador")
        window.geometry("790x650")
        text = tk.Text(window, wrap="word", padx=20, pady=20, font=("Segoe UI", 11))
        scroll = ttk.Scrollbar(window, command=text.yview)
        text.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        text.insert("1.0", HELP)
        text.configure(state="disabled")

    def close(self):
        if self.confirm_discard():
            self.root.destroy()


def parse_view(data):
    view = data.get("camera", {})
    execution = data.get("execucao", {})
    if not isinstance(view, dict) or not isinstance(execution, dict):
        raise ValueError("Configuração de câmera ou execução inválida.")

    def boolean(key, default):
        value = view.get(key, default)
        if not isinstance(value, bool):
            raise ValueError("%s deve ser true ou false." % key)
        return value

    camera = Camera(number(view.get("yaw", -25), "Rotação", -36000, 36000),
                    number(view.get("pitch", 58), "Inclinação da câmera", -90, 90),
                    number(view.get("zoom", 170), "Zoom", 1e-8, 1e12),
                    boolean("perspectiva", True),
                    number(view.get("pan_x", 0), "Deslocamento", -1e8, 1e8),
                    number(view.get("pan_y", 0), "Deslocamento", -1e8, 1e8))
    focus = view.get("foco", "")
    if not isinstance(focus, str):
        raise ValueError("Foco da câmera inválido.")
    return camera, {"size": number(view.get("tamanho", 1), "Tamanho", 1e-4, 1e8),
                    "physical": boolean("tamanho_fisico", False), "focus": focus,
                    "speed": number(execution.get("dias_por_segundo", 10), "Velocidade", 1e-8, 1e8),
                    "step": number(execution.get("passo_max_dias", .05), "Passo", 1e-10, 1e6),
                    "trails": boolean("trilhas", True), "names": boolean("nomes", True),
                    "grid": boolean("grade", True)}


def read_json(path):
    if Path(path).stat().st_size > 2_000_000:
        raise ValueError("Arquivo maior que 2 MB; cenário não suportado.")
    with open(path, encoding="utf-8-sig") as stream:
        return json.load(stream)


def run_tests():
    """Verificações físicas e de persistência independentes de Tk."""
    import unittest
    import copy

    class Tests(unittest.TestCase):
        def binary(self):
            star = Body("Sol", SUN_KG, 696340, [0, 0, 0], [0, 0, 0], "Estrela")
            p, v = orbital_state(star, EARTH_KG, 1, 0, 0, 0, 0, 0)
            planet = Body("Terra", EARTH_KG, 6371, p, v, "Em órbita", star.uid)
            result = System([star, planet], epsilon=0)
            result.center_momentum()
            return result

        def test_circular_orbit_one_period(self):
            s = self.binary()
            initial = sub(s.bodies[1].pos, s.bodies[0].pos)
            period = 2 * math.pi / math.sqrt(G * (SUN_KG + EARTH_KG))
            e0 = s.energy()
            for _ in range(4000):
                s.step(period / 4000)
            self.assertLess(norm(sub(sub(s.bodies[1].pos, s.bodies[0].pos), initial)), 6e-6)
            self.assertLess(abs(s.energy() / e0 - 1), 1e-10)

        def test_momentum_conservation(self):
            s = scenario("Estrelas binárias")
            for _ in range(700):
                s.step(.1)
            total = sum(b.mass for b in s.bodies)
            bulk = [sum(b.mass / total * b.vel[k] for b in s.bodies) for k in range(3)]
            self.assertLess(norm(bulk), 1e-15)
            self.assertLess(norm(s.barycenter()), 1e-12)

        def test_free_motion(self):
            b = Body("Livre", 1000, 0, [1, 2, 3], [.1, -.2, .3])
            s = System([b])
            s.step(2)
            for actual, expected in zip(b.pos, [1.2, 1.6, 3.6]):
                self.assertAlmostEqual(actual, expected)

        def test_time_reversal(self):
            s = self.binary()
            before = copy.deepcopy(s.bodies)
            for _ in range(100):
                s.step(.1)
            for _ in range(100):
                s.step(-.1)
            for a, b in zip(s.bodies, before):
                self.assertLess(norm(sub(a.pos, b.pos)), 1e-12)
                self.assertLess(norm(sub(a.vel, b.vel)), 1e-12)

        def test_inclined_orbit_and_vis_viva(self):
            star = Body("Sol", SUN_KG, 1, [4, 5, 6], [.1, .2, .3], "Estrela")
            p, v = orbital_state(star, EARTH_KG, 2, .3, 90, 0, 0, 90)
            r = norm(sub(p, star.pos))
            self.assertAlmostEqual(p[1], star.pos[1])
            self.assertGreater(abs(p[2] - star.pos[2]), 1)
            self.assertAlmostEqual(dot(sub(v, star.vel), sub(v, star.vel)),
                                   G * (SUN_KG + EARTH_KG) * (2 / r - .5), places=12)

        def test_mutual_attraction(self):
            s = self.binary()
            acc = s.accelerations([b.pos for b in s.bodies])
            self.assertGreater(acc[0][0], 0)
            self.assertLess(acc[1][0], 0)
            self.assertAlmostEqual(abs(acc[0][0] / acc[1][0]), EARTH_KG / SUN_KG)

        def test_collision_is_atomic_and_swept(self):
            a = Body("A", 1, 1, [-.1, 0, 0], [1, 0, 0])
            b = Body("B", 1, 1, [.1, 0, 0], [-1, 0, 0])
            s = System([a, b])
            with self.assertRaises(CollisionError):
                s.step(.2)
            self.assertEqual(a.pos, [-.1, 0, 0])
            self.assertEqual(s.days, 0)

        def test_softening_at_coincidence(self):
            s = System([Body("A", 1, 0, [0, 0, 0], [0, 0, 0]),
                        Body("B", 1, 0, [0, 0, 0], [0, 0, 0])], stop_collision=False)
            s.step(.1)
            self.assertTrue(math.isfinite(s.energy()))

        def test_round_trip(self):
            s = scenario("Órbitas em 3D")
            s.step(.01)
            data = json.loads(json.dumps(system_to_dict(s)))
            other = system_from_dict(data)
            self.assertEqual(other.days, s.days)
            for a, b in zip(s.bodies, other.bodies):
                self.assertEqual(a.uid, b.uid)
                self.assertEqual(a.pos, b.pos)
                self.assertLess(norm(sub(a.vel, b.vel)), 1e-15)

        def test_reject_invalid_scenarios(self):
            data = system_to_dict(self.binary())
            for bad in (float("nan"), float("inf"), -1, "ruim"):
                changed = copy.deepcopy(data)
                changed["corpos"][0]["massa_kg"] = bad
                with self.assertRaises(ValueError):
                    system_from_dict(changed)
            data["corpos"][1]["estrela_id"] = "inexistente"
            with self.assertRaises(ValueError):
                system_from_dict(data)

        def test_visual_scales_leave_dynamics_identical(self):
            s = self.binary()
            altered = copy.deepcopy(s)
            altered.bodies[0].visual_px = 100
            altered.bodies[0].distance_factor = 30
            altered.bodies[1].distance_factor = .02
            altered.bodies[1].visual_px = 200
            for _ in range(10):
                s.step(.1)
                altered.step(.1)
            for a, b in zip(s.bodies, altered.bodies):
                self.assertEqual(a.pos, b.pos)
                self.assertEqual(a.vel, b.vel)

        def test_visual_parent_anchor(self):
            s = self.binary()
            a, b = s.bodies
            a.pos, b.pos = [3, 0, 0], [4, 0, 0]
            a.distance_factor, b.distance_factor = 2, 5
            self.assertEqual(display_position(b, {x.uid: x for x in s.bodies}), [11, 0, 0])

        def test_camera_depth_and_rotation(self):
            c = Camera(yaw=0, pitch=0, zoom=100, perspective=False)
            p = c.project([1, 0, 0], [0, 0, 0], 800, 600)
            self.assertEqual(p[:2], (500, 300))
            c.yaw = 90
            p = c.project([1, 0, 0], [0, 0, 0], 800, 600)
            self.assertAlmostEqual(p[0], 400)
            self.assertAlmostEqual(p[1], 200)
            c.perspective = True
            self.assertIsNone(c.project([0, 0, -100], [0, 0, 0], 800, 600))

        def test_adaptive_step_shrinks_at_approach(self):
            s = self.binary()
            h1 = s.safe_step(1000)
            s.bodies[1].pos = add(s.bodies[0].pos, [.01, 0, 0])
            self.assertLess(s.safe_step(1000), h1 / 100)

        def test_all_examples_evolve(self):
            for name in ("Sistema solar interno", "Órbitas em 3D", "Estrelas binárias", "Corpo livre em passagem", "Vazio"):
                s = scenario(name)
                validate_links(s.bodies)
                for _ in range(30):
                    s.step(.01)
                self.assertTrue(math.isfinite(s.energy()))

        def test_camera_and_run_validation(self):
            c, opts = parse_view({"camera": {"zoom": 200, "tamanho": 5, "perspectiva": False},
                                  "execucao": {"dias_por_segundo": 3, "passo_max_dias": .001}})
            self.assertEqual(c.zoom, 200)
            self.assertFalse(c.perspective)
            self.assertEqual(opts["size"], 5)
            self.assertEqual(opts["step"], .001)
            with self.assertRaises(ValueError):
                parse_view({"camera": {"zoom": -1}})

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Tests)
    return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1


def main():
    parser = argparse.ArgumentParser(description="Simulador Celeste 3D, somente biblioteca padrão.")
    parser.add_argument("--testar", action="store_true", help="executa testes sem abrir a interface")
    parser.add_argument("--cenario", metavar="ARQUIVO", help="abre um cenário JSON (estado físico)")
    parser.add_argument("--gerar-exemplos", metavar="PASTA", help="grava os cenários didáticos em JSON")
    args = parser.parse_args()
    if args.testar:
        return run_tests()
    if args.gerar_exemplos:
        folder = Path(args.gerar_exemplos)
        folder.mkdir(parents=True, exist_ok=True)
        for filename, name in (("sistema_solar", "Sistema solar interno"), ("orbitas_3d", "Órbitas em 3D"),
                               ("binarias", "Estrelas binárias"), ("corpo_livre", "Corpo livre em passagem")):
            (folder / (filename + ".json")).write_text(json.dumps(system_to_dict(scenario(name)),
                                                       ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print("Exemplos gravados em:", folder.resolve())
        return 0
    initial = None
    if args.cenario:
        try:
            initial = system_from_dict(read_json(args.cenario))
        except (OSError, ValueError, TypeError) as exc:
            print("Não foi possível abrir o cenário:", exc, file=sys.stderr)
            return 1
    return launch_gui(initial)


if __name__ == "__main__":
    raise SystemExit(main())
