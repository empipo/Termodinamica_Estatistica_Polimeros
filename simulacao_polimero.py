import copy
import random
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def gerar_polimero(
    rede: np.ndarray,
    tamanho_polimero: int,
    id_polimero: int,
    MAX_TENTATIVAS: int = 10_000,
) -> Optional[list[tuple[int, int]]]:

    altura, largura = rede.shape
    posicoes_livres = list(zip(*np.where(rede == 0)))

    if len(posicoes_livres) == 0:
        return None  # Rede cheia

    for _ in range(MAX_TENTATIVAS):
        y_inicial, x_inicial = random.choice(posicoes_livres)
        polimero = [(x_inicial, y_inicial)]
        rede[y_inicial, x_inicial] = id_polimero
        pos_atual = (x_inicial, y_inicial)

        for _ in range(1, tamanho_polimero):
            vizinhos_possiveis = [
                (pos_atual[0] + 1, pos_atual[1]),
                (pos_atual[0] - 1, pos_atual[1]),
                (pos_atual[0], pos_atual[1] + 1),
                (pos_atual[0], pos_atual[1] - 1),
            ]
            random.shuffle(vizinhos_possiveis)

            passo_valido = False
            for prox_x, prox_y in vizinhos_possiveis:
                if (
                    0 <= prox_x < largura
                    and 0 <= prox_y < altura
                    and rede[prox_y, prox_x] == 0
                ):
                    pos_atual = (prox_x, prox_y)
                    polimero.append(pos_atual)
                    rede[prox_y, prox_x] = id_polimero
                    passo_valido = True
                    break

            if not passo_valido:
                # desfaz
                for px, py in polimero:
                    rede[py, px] = 0
                break

        else:
            return polimero

    return None


def criar_rede(
    TAMANHO_REDE_X,
    TAMANHO_REDE_Y,
    NUMERO_DE_POLIMEROS,
    TAMANHO_DO_POLIMERO,
    MAX_TENTATIVAS,
):
    for i in range(MAX_TENTATIVAS):
        rede = np.zeros((TAMANHO_REDE_Y, TAMANHO_REDE_X))
        polimeros = {}
        id_atual = 1

        while id_atual <= NUMERO_DE_POLIMEROS:
            coords = gerar_polimero(rede, TAMANHO_DO_POLIMERO, id_atual, MAX_TENTATIVAS)
            if coords is None:
                break
            polimeros[id_atual] = coords
            id_atual += 1

        if id_atual > NUMERO_DE_POLIMEROS:
            return rede, polimeros

    raise RuntimeError(
        "Não foi possível gerar uma rede completa após muitas tentativas."
    )


def calcula_energia_sistema(
    rede: list[list],
    polimeros_gerados: dict[int, list[list]],
    EAl=1,
    EAd=1,
    EBl=1,
    EBd=1,
):
    Energia_total = 0
    dicionario = dict()
    lista_energia = []

    for chave, valor in polimeros_gerados.items():
        passo = 0
        dicionario[chave] = {
            "frontais": 0,
            "diagonais": 0,
            "diferentes frontais": 0,
            "diferentes diagonais": 0,
        }
        energia_iteracao = 0

        for v in valor:
            rede_copia = np.copy(rede)
            y, x = v
            rede_copia[v[1], v[0]] = 0

            if passo != 0:
                ant = valor[passo - 1]
                rede_copia[ant[1], ant[0]] = 0

            if passo != len(valor) - 1:
                prox = valor[passo + 1]
                rede_copia[prox[1], prox[0]] = 0

            passo += 1

            rede_pad = np.pad(
                rede_copia, pad_width=1, mode="constant", constant_values=0
            )
            campo = rede_pad[x : x + 3, y : y + 3]

            # Interações Laterais
            vizinhos_possiveis_laterais = [(0, 1), (1, 0), (1, 2), (2, 1)]
            for coord in vizinhos_possiveis_laterais:
                campo_valor = campo[coord]
                if campo_valor == 0:
                    pass
                elif campo_valor == chave:
                    Energia_total += EAl
                    energia_iteracao += EAl
                    dicionario[chave]["frontais"] += 1
                else:
                    Energia_total += EBl
                    energia_iteracao += EBl
                    dicionario[chave]["diferentes frontais"] += 1

            # Interações Diagonais
            vizinhos_possiveis_diagonais = [(0, 0), (0, 2), (2, 0), (2, 2)]
            for coord_d in vizinhos_possiveis_diagonais:
                campo_valor = campo[coord_d]
                if campo_valor == 0:
                    pass
                elif campo_valor == chave:
                    Energia_total += EAd
                    energia_iteracao += EAd
                    dicionario[chave]["diagonais"] += 1
                else:
                    Energia_total += EBd
                    energia_iteracao += EBd
                    dicionario[chave]["diferentes diagonais"] += 1

        dicionario[chave]["frontais"] /= 2
        dicionario[chave]["diagonais"] /= 2
        lista_energia.append(energia_iteracao / 2)

    Energia_total = Energia_total / 2
    return dicionario, rede, Energia_total, lista_energia


def passo_metro_polimero(
    rede: list[list],
    polimeros_gerados: dict[int, list[list]],
    T: int,
    TAMANHO_DO_POLIMERO: int,
    MAX_TENTATIVAS: int = 10_000,
    EAl=1,
    EAd=1,
    EBl=1,
    EBd=1,
):
    polimeros_originais = copy.deepcopy(polimeros_gerados)
    _, _, E_antiga, _ = calcula_energia_sistema(
        rede, polimeros_originais, EAl, EAd, EBl, EBd
    )

    id_escolhido = random.choice(list(polimeros_gerados.keys()))
    polimero_antigo = polimeros_gerados[id_escolhido]

    rede_nova = copy.deepcopy(rede)
    for i in range(len(rede_nova)):
        for j in range(len(rede_nova[i])):
            if rede_nova[i][j] == id_escolhido:
                rede_nova[i][j] = 0

    coordenadas = gerar_polimero(
        rede_nova, TAMANHO_DO_POLIMERO, id_escolhido, MAX_TENTATIVAS
    )

    if coordenadas:
        polimeros_gerados[id_escolhido] = coordenadas
        for x, y in coordenadas:
            rede_nova[y][x] = id_escolhido
    else:
        return rede, polimeros_originais

    _, _, E_nova, _ = calcula_energia_sistema(
        rede_nova, polimeros_gerados, EAl, EAd, EBl, EBd
    )
    dE = E_nova - E_antiga

    probabilidade = np.random.rand()
    if dE <= 0 or probabilidade < np.exp(-dE / T):
        return rede_nova, polimeros_gerados
    else:
        polimeros_gerados[id_escolhido] = polimero_antigo
        return rede, polimeros_originais


def calcula_distancia_ponta_a_ponta(polimeros_gerados):
    distancias_r = {}
    for id_polimero, coordenadas in polimeros_gerados.items():
        if len(coordenadas) < 2:
            distancias_r[id_polimero] = 0.0
            continue
        x_inicio, y_inicio = coordenadas[0]
        x_fim, y_fim = coordenadas[-1]
        distancia_r = np.sqrt((x_fim - x_inicio) ** 2 + (y_fim - y_inicio) ** 2)
        distancias_r[id_polimero] = distancia_r
    return distancias_r


def calcula_capacidade_calorifica(energias: list[float], T: float) -> float:
    if T == 0:
        return 0.0
    energias_np = np.array(energias)
    media_E_quadrado = np.mean(energias_np**2)
    quadrado_media_E = np.mean(energias_np) ** 2
    flutuacao_energia = media_E_quadrado - quadrado_media_E
    Cv = flutuacao_energia / (T**2)
    return Cv


def simular(
    TAMANHO_REDE_X,
    TAMANHO_REDE_Y,
    NUMERO_DE_POLIMEROS,
    TAMANHO_DO_POLIMERO,
    temps,
    equil_steps,
    prod_steps,
    EAl,
    EAd,
    EBl,
    EBd,
    MAX_TENTATIVAS=10_000,
):

    Energia_media = []
    Distancia_R_media = []
    Cv_media = []

    for T in temps:
        rede, polimeros_gerados = criar_rede(
            TAMANHO_REDE_X,
            TAMANHO_REDE_Y,
            NUMERO_DE_POLIMEROS,
            TAMANHO_DO_POLIMERO,
            MAX_TENTATIVAS,
        )

        for _ in range(equil_steps):
            rede, polimeros_gerados = passo_metro_polimero(
                rede,
                polimeros_gerados,
                T,
                TAMANHO_DO_POLIMERO,
                MAX_TENTATIVAS,
                EAl,
                EAd,
                EBl,
                EBd,
            )

        energias = []
        distancias_R_passo = []

        for _ in range(prod_steps):
            rede, polimeros_gerados = passo_metro_polimero(
                rede,
                polimeros_gerados,
                T,
                TAMANHO_DO_POLIMERO,
                MAX_TENTATIVAS,
                EAl,
                EAd,
                EBl,
                EBd,
            )

            _, _, e, _ = calcula_energia_sistema(
                rede, polimeros_gerados, EAl, EAd, EBl, EBd
            )
            energias.append(e)

            distancias_r_polimeros = calcula_distancia_ponta_a_ponta(polimeros_gerados)
            media_r_passo = np.mean(list(distancias_r_polimeros.values()))
            distancias_R_passo.append(media_r_passo)

        Energia_media.append(np.mean(energias))
        Distancia_R_media.append(np.mean(distancias_R_passo))
        Cv_media.append(calcula_capacidade_calorifica(energias, T))

    return Energia_media, Distancia_R_media, Cv_media


if __name__ == "__main__":
    T_foco = np.linspace(1, 50, 15)
    T_alta = np.linspace(50, 100, 5)
    temps_novo = np.unique(np.concatenate([T_foco, T_alta]))

    equil_steps = 5000
    prod_steps = 1000

    cenarios = {
        "Padrão (Colapso + Segregação)": {"EAl": -1, "EAd": -1, "EBl": 1, "EBd": 1},
        "Sem Energia Intracadeia": {"EAl": 0, "EAd": 0, "EBl": 1, "EBd": 1},
        "Sem Repulsão Intercadeia": {"EAl": -1, "EAd": -1, "EBl": 0, "EBd": 0},
        "Apenas Volume Excluído (SAW)": {"EAl": 0, "EAd": 0, "EBl": 0, "EBd": 0},
    }

    resultados = {}

    print(f"Iniciando simulação comparativa para {len(cenarios)} cenários...")

    for nome_cenario, params in tqdm(cenarios.items(), desc="Cenários"):
        E, R, Cv = simular(
            TAMANHO_REDE_X=12,
            TAMANHO_REDE_Y=12,
            NUMERO_DE_POLIMEROS=10,
            TAMANHO_DO_POLIMERO=10,
            temps=temps_novo,
            equil_steps=equil_steps,
            prod_steps=prod_steps,
            MAX_TENTATIVAS=10_000,
            **params,
        )
        resultados[nome_cenario] = (E, R, Cv)

    fig, axs = plt.subplots(3, 1, figsize=(10, 14), sharex=True)
    fig.suptitle("Comparação de Parâmetros Energéticos em Polímeros", fontsize=16)

    estilos = ["o-", "s-", "^-", "D-"]

    for i, (nome, (E, R, Cv)) in enumerate(resultados.items()):
        estilo = estilos[i % len(estilos)]

        axs[0].plot(temps_novo, E, estilo, markersize=5, label=nome, alpha=0.8)

        axs[1].plot(temps_novo, R, estilo, markersize=5, label=nome, alpha=0.8)

        axs[2].plot(temps_novo, Cv, estilo, markersize=5, label=nome, alpha=0.8)

    axs[0].set_ylabel(r"$\langle E \rangle$ (Energia Média)", fontsize=12)
    axs[0].set_title(r"Energia vs Temperatura", fontsize=14)
    axs[0].grid(True, linestyle="--", alpha=0.5)
    axs[0].legend(fontsize=10)

    axs[1].set_ylabel(r"$\langle R \rangle$ (Distância Ponta a Ponta)", fontsize=12)
    axs[1].set_title(r"Conformação ($\langle R \rangle$) vs Temperatura", fontsize=14)
    axs[1].grid(True, linestyle="--", alpha=0.5)

    axs[2].set_xlabel(r"$T$ (Temperatura)", fontsize=12)
    axs[2].set_ylabel(r"$C_v$ (Capacidade Calorífica)", fontsize=12)
    axs[2].set_title(r"Capacidade Calorífica ($C_v$) vs Temperatura", fontsize=14)
    axs[2].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("comparacao_energias.png", dpi=300)
    print("\nGráfico salvo como 'comparacao_energias.png'")
    plt.show()
