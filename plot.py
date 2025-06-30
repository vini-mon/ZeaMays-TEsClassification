import matplotlib.pyplot as plt
import numpy as np

# --- Dados ---
elementos = ['LTR', 'TIR', 'Helitron', 'MITE', 'LINE', 'SINE']
contagens = [38795, 18995, 16337, 15126, 13046, 5272]

# --- Criação do Gráfico ---

# Define o tamanho da figura para melhor visualização
plt.figure(figsize=(10, 6))

# Cria as barras
barras = plt.bar(elementos, contagens, color=['#8B0000', '#4682B4', '#2E8B57', '#DAA520', '#4B0082', '#D2691E'])

# --- Personalização ---

# Adiciona título e rótulos aos eixos
plt.title('Distribuição de Elementos Transponíveis', fontsize=16)
plt.xlabel('Elementos Transponíveis', fontsize=12)
plt.ylabel('Contagem', fontsize=12)

# Adiciona os valores no topo de cada barra para clareza
for barra in barras:
    yval = barra.get_height()
    plt.text(barra.get_x() + barra.get_width()/2.0, yval, int(yval), va='bottom', ha='center', fontsize=10) # va='bottom' para colocar acima

# Garante que tudo se ajuste bem na imagem final
plt.tight_layout()

# --- Exibição do Gráfico ---
plt.show()