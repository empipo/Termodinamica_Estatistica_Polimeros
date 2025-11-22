<h1 align="center"> Mecânica Estatística para entender a formação de cadeias poliméricas </h1>
<h2 align="center"> Ensemble Canônico e Simulação</h2> 

<p align="center">
   Polémeros 
  &nbsp;&bull;&nbsp; Ensemble Canônico
  &nbsp;&bull;&nbsp; Método do Monte Carlo
  &nbsp;&bull;&nbsp; Andar do Bêbado
</p>

Autores: Ana Luz Pereira Mendes, Emanuel Piveta Pozzobon , Pedro Coelho G. de Freitas 

Orientação: Prof. Dr. F. Crasto de Lima

-----------
<p align="center">
<img loading="lazy" src="http://img.shields.io/static/v1?label=STATUS&message=EM%20DESENVOLVIMENTO&color=GREEN&style=for-the-badge"/>
</p>

Este repositório contém a implementação de um simulador estocástico para estudar a configuração topológica e o comportamento energético de sistemas poliméricos em uma rede bidimensional. O projeto utiliza o ensemble canônico, o algoritmo do Andar do Bêbado e o Método de Monte Carlo para explorar conformações possíveis, avaliar energias e analisar a influência da temperatura na organização espacial dos polímeros.

# Objetivo 📌 

Investigar, por meio de simulações computacionais, o comportamento de N homopolímeros distribuídos em uma rede 2D. O simulador busca:

- Determinar conformações poliméricas possíveis em uma rede usando caminhadas aleatórias.

- Calcular a energia total do sistema com base em interações entre primeiros vizinhos.

- Utilizar o método de Monte Carlo para aceitar/rejeitar novas configurações conforme a distribuição do ensemble canônico.

- Analisar como a temperatura afeta a probabilidade de microestados, a energia do sistema e o enovelamento das cadeias.

# Contexto Teórico

Polímeros podem ser modelados como macromoléculas compostas por unidades repetitivas ligadas covalentemente. Sua topologia e dinâmica configuracional estão presentes em fenômenos biológicos (como conformação de DNA e parede celular), bem como em materiais sintéticos (biodegradáveis ou condutores).

Para sistemas cujo estado energético varia dinamicamente, o ensemble canônico é uma poderosa ferramenta estatística. A probabilidade de um sistema estar em um microestado j, de energia 𝐸_j, é dada por:

$$
P_j = \frac{e^{-\beta E_j}}{Z}
$$
com
$$
\beta = \frac{1}{k_B T}
$$

A função de partição 𝑍 é essencial para normalizar esta distribuição.

Novas configurações são geradas neste projeto pelo Andar do Bêbado, uma caminhada aleatória simples que modela a distribuição espacial dos monômeros.
A simulação utiliza o Método de Monte Carlo para aceitar alterações estruturais com base no peso de Boltzmann.

---------------------



