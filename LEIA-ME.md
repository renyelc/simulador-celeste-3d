# Simulador Celeste 3D

Aplicativo desktop em **Python 3.9 ou superior**, somente com a biblioteca padrão. Não usa NumPy, Matplotlib, Pygame, OpenGL ou instalação por `pip`. Funciona offline.

## Começar no Windows

1. Extraia o ZIP para uma pasta.
2. Abra `Iniciar_Windows.bat` com um duplo clique.
3. Clique em **Iniciar**. O sistema solar interno já estará carregado.

Também pode executar no terminal, dentro da pasta:

```console
python simulador_celeste.py
```

Se o comando do seu Windows for `py`, use:

```console
py -3 simulador_celeste.py
```

O arquivo `simulador_celeste.py` é autossuficiente: os exemplos e a ajuda também estão embutidos nele. Pode copiá-lo sozinho para outra pasta.

**Requisito gráfico:** Python com Tkinter e Tcl/Tk, normalmente incluídos na instalação oficial para Windows. Algumas distribuições mínimas de Python em Linux/macOS não incluem Tk. Nesse caso, a interface não abre sem um Python que já tenha esse componente. Os testes numéricos funcionam sem Tk e sem tela. Para conferir seu Python, execute `python -m tkinter`: deve abrir uma janela de teste.

## O que o aplicativo faz

| Recurso | Funcionamento |
|---|---|
| Múltiplos corpos | Todos exercem gravidade sobre todos; as estrelas também se movem. |
| Movimento em três dimensões | Cada corpo tem posição X/Y/Z e velocidade Vx/Vy/Vz. |
| Visualização 3D | Projeção em Canvas com rotação, perspectiva, profundidade e discos sombreados. |
| Tipos de corpo | Livre, Estrela e Em órbita. |
| Órbitas iniciais | Gerador de órbitas circulares ou elípticas, inclusive inclinadas. |
| Escalas independentes | Distância em pixels/UA e tamanho visual ajustados separadamente. |
| Ajuste por corpo | Raio visual em pixels e fator visual da distância. |
| Controle temporal | Iniciar, pausar, avançar um passo, reiniciar e controlar dias/s. |
| Inspeção | Trilhas, eixos, grade, foco da câmera, distâncias físicas e energia. |
| Cenários | Salvar e abrir arquivos JSON com o estado físico e os ajustes visuais. |
| Colisões | Pausa antes de um passo com contato, usando os raios físicos. |

## Como controlar a câmera

- **Botão esquerdo + arrastar:** girar a vista.
- **Roda do mouse:** aumentar ou reduzir a escala das distâncias.
- **Botão direito + arrastar**, ou **Shift + esquerdo:** deslocar a vista.
- **Clique no corpo:** selecionar; também pode selecionar na lista.
- **Duplo clique na lista:** editar o corpo.
- **XY, XZ e YZ:** observar os planos coordenados. **3D:** voltar à vista oblíqua.
- **Enquadrar:** ajustar as distâncias para o sistema caber na vista, em torno do foco atual.
- **Centro da câmera:** origem, baricentro físico ou um corpo específico.

O painel lateral tem uma barra de rolagem. Use-a para alcançar os controles inferiores em telas menores.

## Adicionar uma estrela e um planeta

1. Selecione o cenário **Vazio** e clique em **Carregar exemplo**.
2. Clique em **+ Novo** e cadastre:
   - Nome: `Minha estrela`.
   - Tipo: **Estrela**.
   - Massa: `1.98847e30` kg.
   - Raio físico: `696340` km.
   - Raio visual: `15` px.
   - Posição X/Y/Z: `0`, `0`, `0` UA.
   - Velocidade X/Y/Z: `0`, `0`, `0` km/s.
3. Salve a estrela e clique novamente em **+ Novo**.
4. Cadastre `Meu planeta`, do tipo **Em órbita**, e selecione a estrela criada.
5. Use massa `5.9722e24` kg, raio físico `6371` km e raio visual `7` px.
6. Abra a aba **Gerador de órbita**. Use:
   - Semieixo maior: `1` UA.
   - Excentricidade: `0` para circular, ou `0.2` para elíptica.
   - Inclinação: `35` graus para enxergar a órbita fora do plano XY.
   - Nó ascendente, argumento do periastro e anomalia verdadeira: `0`.
7. Clique em **Gerar posição e velocidade**. Confira os campos preenchidos na outra aba e clique em **Salvar corpo**.
8. Clique em **Referencial do baricentro** se desejar retirar o movimento global do conjunto. Isso é uma mudança de referencial: preserva posições e velocidades relativas.
9. Clique em **Enquadrar** e depois em **Iniciar**.

**Atenção:** preencher a aba do gerador não basta. É preciso clicar em **Gerar posição e velocidade** e depois em **Salvar corpo**. Ao editar um corpo já existente, os campos do gerador começam com valores padrão; eles não são uma estimativa dos elementos da órbita atual.

O gerador soma a posição e a velocidade da estrela. Ele utiliza a massa da estrela **e** do corpo para calcular a velocidade relativa. A estrela não é reposicionada automaticamente nessa operação; o botão do baricentro pode ser usado depois.

## Corpos livres, estrelas e corpos em órbita

**Livre:** você define a posição e a velocidade. Pode representar um cometa em passagem, um asteroide ou até um planeta em torno de um sistema binário. “Livre” não significa livre de gravidade.

**Estrela:** também recebe posição e velocidade, sofre a gravidade dos demais corpos e pode servir de referência ao gerador de órbitas. Não fica fixa no centro.

**Em órbita:** possui uma estrela associada para o gerador, a medida de distância e o ajuste visual. A trajetória não é imposta: o corpo pode ser perturbado, colidir ou escapar. Também é possível digitar suas coordenadas manualmente.

A associação é sempre a uma estrela. Para simular uma lua em torno de um planeta, cadastre-a como **Livre** e informe manualmente a posição e a velocidade globais adequadas. Este gerador não cria automaticamente órbitas ao redor de planetas.

Se excluir uma estrela, seus corpos associados passam a ser livres. Para mudar uma estrela para outro tipo, primeiro reassocie ou torne livres os corpos que dependem dela.

## Escalas visuais: tamanho e distância separados

### Distância global

O campo **Distância: pixels / UA** controla o espaçamento na tela. Digite o valor e clique em **Aplicar**. A roda do mouse altera esse mesmo valor.

### Tamanho global

O campo **Tamanho global** é um multiplicador independente. Digite, por exemplo, `2` e clique em **Aplicar** para dobrar os raios visuais sem alterar as posições.

- **Manual:** utiliza o raio visual em pixels de cada corpo, multiplicado pelo tamanho global. É o padrão didático, para que planetas sejam visíveis.
- **Físico:** converte o raio físico em km para a mesma escala usada pelas posições e aplica o multiplicador de tamanho. Com multiplicador `1` e fatores de distância `1`, aproxima as proporções físicas. Existe um marcador mínimo de 1,5 px para corpos muito pequenos continuarem selecionáveis. A perspectiva ainda altera o tamanho aparente com a profundidade.

O desenho limita raios extremamente grandes para proteger a renderização; esse caso não representa uma esfera fielmente em escala.

### Ajustes por corpo

Abra **Editar**:

- **Raio visual manual (px):** tamanho individual no modo Manual.
- **Fator distância visual:** `1` preserva a posição visual, `2` duplica seu afastamento e `0.5` o reduz à metade.

Para **Em órbita**, esse afastamento é medido em relação à posição visual da estrela associada. Para estrelas e corpos livres, é medido a partir da origem. Se mudar o fator da estrela, os corpos associados acompanham sua nova posição visual, com seus próprios afastamentos.

Exemplo: um planeta a `0.1 UA` da estrela pode receber fator `10` e aparecer a uma distância visual equivalente a `1 UA`. A gravidade e a distância exibida no painel continuam usando `0.1 UA`.

**Nenhuma escala visual altera massa, raio físico, posição, velocidade, força ou período orbital.** O painel continua mostrando distâncias físicas, de centro a centro. A grade indica coordenadas físicas e não mede diretamente afastamentos distorcidos por fatores individuais. O baricentro físico também não é recalculado usando posições visuais.

As trilhas representam posições globais passadas. Ao seguir uma estrela em movimento, a trilha pode não parecer uma elipse centrada nela. Editar corpos ou mudar o referencial limpa as trilhas.

## Unidades e condições mecânicas

| Dado | Unidade | Exemplo |
|---|---|---|
| Massa | kg | Terra: `5.9722e24` |
| Raio físico | km | Terra: `6371` |
| Posição X/Y/Z | UA | Terra aproximadamente a `1` UA do Sol |
| Velocidade Vx/Vy/Vz | km/s | Aproximadamente `29.78` para órbita circular a 1 UA |
| Tempo e passo | dias | `0.05` dia corresponde a 72 minutos |
| Inclinações e ângulos | graus | `90` gera um plano perpendicular ao XY |

Uma UA equivale a **149.597.870,7 km**. Aceita ponto ou vírgula decimal e notação científica, como `1,98847e30`. Não use separadores de milhares.

No gerador:

- **Semieixo maior a:** escala da elipse; numa órbita circular é o raio orbital.
- **Excentricidade e:** de `0` até `0.99`; `0` é circular. Periastro = `a(1−e)`; apoastro = `a(1+e)`.
- **Inclinação i:** ângulo do plano orbital em relação ao plano XY.
- **Longitude do nó ascendente:** orientação da linha de interseção entre o plano orbital e o plano XY.
- **Argumento do periastro:** orientação da menor distância dentro do plano orbital.
- **Anomalia verdadeira:** posição angular inicial medida a partir do periastro.

## Integração e precisão

O programa integra a gravitação newtoniana com **velocity Verlet**. Para cada par, a aceleração usa:

```text
a_i = G · m_j · (r_j − r_i) / (|r_j − r_i|² + ε²)^(3/2)
```

As posições internas estão em UA, o tempo em dias e as massas em kg. A constante gravitacional é convertida a essas unidades a partir de `G = 6.67430e−11 m³/(kg·s²)`.

O passo escolhido é limitado pelo máximo informado e por 3% das escalas locais de queda livre e cruzamento dos pares. Aproximações exigem passos menores. Velocity Verlet com passo fixo é simplético; a adaptação do passo perde essa propriedade estrita.

**Dias/s é uma meta de reprodução.** O programa limita o trabalho por ciclo para manter a interface utilizável. Quando o processamento fica pesado, avança menos tempo físico e mostra a velocidade efetiva na barra inferior.

Comece com `10` dias/s e passo máximo de `0.05` dia. Para encontros próximos, órbitas muito pequenas ou objetos muito massivos, reduza o passo. Use o desvio de energia como um diagnóstico; compare também o resultado ao reduzir o passo pela metade. A estabilidade aparente de uma animação, sozinha, não garante precisão.

**Suavização:** o padrão é `1e−7 UA`, aproximadamente 15 km. Evita singularidades e modifica a gravitação em distâncias muito pequenas. Para escalas pequenas, ajuste esse valor ao problema. Valor zero usa a lei de Newton sem suavização, mas corpos pontuais coincidentes tornam a força indefinida.

**Colisões:** o programa pausa antes de um passo que teria contato entre raios físicos. A interseção é aproximada por segmentos percorridos durante o passo. Não modela fusão, fragmentação, deformação ou ricochete. Se desmarcar a pausa, os corpos podem se atravessar sob a gravitação suavizada.

## Salvar, abrir e reiniciar

- **Salvar:** grava o estado físico atual, o tempo decorrido, a câmera, os ajustes de exibição e de execução em JSON. Não grava o histórico de trilhas.
- **Abrir:** restaura esse estado pausado e o usa como nova referência para Reiniciar.
- **Reiniciar:** volta à última referência estabelecida; não necessariamente ao exemplo original.
- Adicionar, editar, excluir, alterar suavização/colisões ou mudar o referencial cria uma nova referência no instante atual.
- Para recuperar o exemplo original, use **Carregar exemplo**.
- O JSON é lido apenas como dados; não é executado como código.

Também pode abrir um arquivo pelo terminal:

```console
python simulador_celeste.py --cenario exemplos/orbitas_3d.json
```

Essa opção de terminal importa o estado físico e enquadra a cena. Para recuperar também os ajustes de câmera de um arquivo salvo, use **Abrir** dentro do aplicativo.

## Exemplos incluídos

1. **Sistema solar interno:** Sol, Mercúrio, Vênus, Terra e Marte. Elementos orbitais aproximados e fases ilustrativas; não corresponde a uma data real.
2. **Órbitas em 3D:** três planetas com inclinações de 0°, 45° e 80° e diferentes excentricidades.
3. **Estrelas binárias:** duas estrelas que orbitam o baricentro e um planeta circumbinário com condição inicial aproximada. O planeta é cadastrado como Livre e sofre a força de ambas.
4. **Corpo livre em passagem:** Sol, Terra e um visitante que sofre deflexão gravitacional.

Para recriar os arquivos JSON:

```console
python simulador_celeste.py --gerar-exemplos exemplos
```

## Verificação e limites

Execute os 16 testes incluídos:

```console
python simulador_celeste.py --testar
```

Cobrem órbita circular completa, energia, conservação do momento e baricentro, reversibilidade temporal, movimento livre, atração mútua, órbita inclinada, relação vis-viva, colisão, suavização, JSON, validação de entradas, projeção 3D, passo adaptativo e independência das escalas visuais. Veja `VALIDACAO.txt` para o resultado obtido durante a preparação.

**Validação realizada:** os testes numéricos e de persistência passaram em Python 3.12. A sintaxe foi compilada e os quatro arquivos de exemplo foram carregados e evoluídos. O ambiente de preparação não possui sessão gráfica; não foi possível abrir a janela de Tkinter nem verificar visualmente seu layout. Os gestos do mouse e diálogos devem ser conferidos na primeira execução no seu computador.

É um simulador didático newtoniano, sem relatividade, marés, rotação axial, atmosfera, pressão, evolução estelar ou efemérides astronômicas. A visualização usa projeção matemática em Canvas, com discos sombreados e ordenação aproximada por profundidade; não é um renderizador volumétrico de esferas.

As forças têm custo O(N²). Comece com 5 a 30 corpos; o desempenho depende do computador. O limite de cadastro é 150 corpos e o arquivo de cenário deve ter até 2 MB. Não há dependências de terceiros nem acesso à rede durante o uso.

## Arquivos do pacote

```text
simulador_celeste.py    Programa completo e testes integrados
Iniciar_Windows.bat    Atalho de execução para Windows
LEIA-ME.md             Este guia
VALIDACAO.txt          Resultado dos testes e limites da verificação
exemplos/              Quatro cenários JSON
```

O código está organizado em física (`Body`, `System`, `orbital_state`), visualização (`Camera`, `display_position`), interface (`Application`, `BodyDialog`) e validação (`run_tests`). Todos estão no mesmo arquivo para facilitar transporte e execução.
