# Simulador Celeste 3D — versão 1.0

Simulador de dinâmicas gravitacionais de corpos celestes.

Aplicativo desktop em Python 3.9+ com Tkinter, usando apenas a biblioteca padrão, sem instalação por pip.

## Executar

Baixe o repositório em **Code → Download ZIP**, extraia os arquivos e abra `Iniciar_Windows.bat` no Windows. Ou execute:

```bash
python simulador_celeste.py
```

Requer uma instalação de Python com Tkinter/Tcl-Tk. Para conferir: `python -m tkinter`.

## Recursos

- Gravitação newtoniana entre todos os corpos, com posições e velocidades em três dimensões.
- Cadastro de corpos livres, estrelas e corpos associados a estrelas.
- Geração de condições iniciais para órbitas circulares ou elípticas e inclinadas.
- Projeção 3D com câmera interativa, trilhas e escalas visuais independentes de distância e tamanho.
- Integração temporal velocity Verlet, redução heurística do passo e pausa em colisões.
- Salvamento e abertura de cenários JSON; quatro exemplos incluídos.

## Documentação

- [Guia completo de uso](LEIA-ME.md)
- [Explicação do código e do motor temporal — Word](docs/Simulador_Celeste_3D_Documentacao_v1_0.docx)
- [Resultados e limites da validação](VALIDACAO.txt)
- [Histórico de versões](CHANGELOG.md)

## Testes

```bash
python simulador_celeste.py --testar
```

São 16 testes integrados, executáveis sem abrir a interface. Todos passaram novamente na preparação desta publicação. A janela gráfica não foi verificada visualmente no ambiente de preparação.

## Limites

Simulador didático, sem relatividade ou efemérides astronômicas. As forças têm custo O(N²). A visualização usa projeção em Canvas e discos sombreados. As escalas visuais não modificam a dinâmica física. Consulte o guia para os limites do passo adaptativo e das colisões.
