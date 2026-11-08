# -*- coding: utf-8 -*-
"""
Gera as figuras didáticas (PNG) usadas na documentação de séries temporais.
Ambiente: Python 3.12 + numpy + matplotlib (sem internet).
As simulações têm sementes fixas para reprodutibilidade.
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({
    "figure.dpi": 130,
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

OUT = "/workspace/st_doc/figuras"
os.makedirs(OUT, exist_ok=True)

AZUL = "#1f4e79"
LARANJA = "#d55e00"
VERDE = "#0072b2"
CINZA = "#666666"


def fmt_br(x, pos=None):
    return f"{x:.2f}".replace(".", ",")


# ----------------------------------------------------------------------
# FIGURA 1 — Componentes de uma série temporal (decomposição clássica)
# ----------------------------------------------------------------------
def figura_componentes():
    n = 60  # 5 anos mensais
    t = np.arange(1, n + 1)
    rng = np.random.default_rng(42)

    tendencia = 100 + 1.1 * t - 0.02 * (t - 30) ** 2 / 6  # leve curvatura
    sazonal = 9 * np.sin(2 * np.pi * t / 12) + 3 * np.sin(2 * np.pi * t / 4)
    ruido = rng.normal(0, 2.2, size=n)
    y = tendencia + sazonal + ruido

    # Decomposição clássica (aditiva) — passo a passo
    m = 12
    # Passo 1: média móvel centrada 2x12 para estimar tendência-ciclo
    def ma(v, w):
        pad = np.pad(v, (w, w), mode="edge")
        return np.convolve(pad, np.ones(2 * w) / (2 * w), mode="same")[w:-w]
    T_est = ma(y, m)  # 2x12

    # Passo 2: série destendenciada
    det = y - T_est

    # Passo 3: média por estação (meses) e ajuste para soma zero
    S_est = np.zeros(m)
    counts = np.zeros(m)
    for i in range(n):
        S_est[i % m] += det[i]
        counts[i % m] += 1
    S_est /= counts
    S_est -= S_est.mean()
    S_full = np.array([S_est[i % m] for i in range(n)])

    # Passo 4: restante
    R_est = y - T_est - S_full

    fig, axs = plt.subplots(4, 1, figsize=(9, 8.6), sharex=True)
    axs[0].plot(t, y, color=AZUL, lw=1.6, label="Série observada (Y)")
    axs[0].set_title("① Série observada — tendência + sazonalidade + ruído")
    axs[0].legend(loc="upper left", frameon=False)
    axs[0].set_ylabel("valor")

    axs[1].plot(t, T_est, color=LARANJA, lw=2.2, label="Tendência-ciclo estimada")
    axs[1].set_title("② Tendência-ciclo (média móvel centrada 2×12)")
    axs[1].legend(loc="upper left", frameon=False)
    axs[1].set_ylabel("valor")

    axs[2].plot(t, S_full, color=VERDE, lw=1.6, label="Componente sazonal")
    axs[2].set_title("③ Componente sazonal (padrão repetido mês a mês)")
    axs[2].legend(loc="upper left", frameon=False)
    axs[2].set_ylabel("desvio")

    axs[3].plot(t, R_est, color=CINZA, lw=1.2, label="Restante (ruído)")
    axs[3].set_title("④ Restante / irregular (o que sobrou após retirar ② e ③)")
    axs[3].legend(loc="upper left", frameon=False)
    axs[3].set_xlabel("tempo (meses)")
    axs[3].set_ylabel("desvio")

    fig.suptitle("Decomposição clássica de uma série temporal simulada (exemplo didático)",
                 fontsize=12.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(f"{OUT}/fig1_componentes.png", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# FIGURA 2 — Assinaturas de ACF e PACF: AR(1) vs MA(1)
# ----------------------------------------------------------------------
def acf_empirica(x, lag_max):
    x = x - x.mean()
    s = np.sum(x ** 2)
    out = []
    for k in range(1, lag_max + 1):
        out.append(np.sum(x[: len(x) - k] * x[k:]) / s)
    return np.array(out)


def pacf_empirica(x, lag_max):
    """PACF via OLS: coeficiente de y_{t-k} na regressão de y_t em y_{t-1}..y_{t-k}."""
    x = x - x.mean()
    out = []
    for k in range(1, lag_max + 1):
        X = np.column_stack([x[k - j - 1 : len(x) - j - 1] for j in range(k)])
        y = x[k:]
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        out.append(beta[-1])
    return np.array(out)


def figura_acf_pacf():
    rng = np.random.default_rng(7)
    n = 800

    # AR(1): y_t = 0.8*y_{t-1} + e_t
    e = rng.normal(0, 1, n + 100)
    ar = np.zeros(n + 100)
    for i in range(1, n + 100):
        ar[i] = 0.8 * ar[i - 1] + e[i]
    ar = ar[100:]

    # MA(1): y_t = e_t + 0.8*e_{t-1}
    ma = e[1:] + 0.8 * e[:-1]
    ma = ma[:n]

    lag = 12
    acf_ar = acf_empirica(ar, lag)
    pacf_ar = pacf_empirica(ar, lag)
    acf_ma = acf_empirica(ma, lag)
    pacf_ma = pacf_empirica(ma, lag)

    bands = 1.96 / np.sqrt(n)

    fig, axs = plt.subplots(2, 2, figsize=(9.5, 6.4), sharex=True)
    lags = np.arange(1, lag + 1)

    def barras(ax, vals, cor, titulo):
        ax.bar(lags, vals, color=cor, width=0.6)
        ax.axhline(bands, color="#999999", ls="--", lw=1)
        ax.axhline(-bands, color="#999999", ls="--", lw=1)
        ax.set_title(titulo)
        ax.set_xticks(lags)
        ax.set_ylim(-1, 1)

    barras(axs[0, 0], acf_ar, AZUL, "AR(1): ACF se extingue aos poucos (cauda)")
    barras(axs[0, 1], pacf_ar, LARANJA, "AR(1): PACF corta após o lag 1 (picos)")
    barras(axs[1, 0], acf_ma, VERDE, "MA(1): ACF corta após o lag 1 (picos)")
    barras(axs[1, 1], pacf_ma, CINZA, "MA(1): PACF se extingue aos poucos (cauda)")

    axs[0, 0].set_ylabel("correlação")
    axs[1, 0].set_ylabel("correlação")
    axs[1, 0].set_xlabel("defasagem (lags)")
    axs[1, 1].set_xlabel("defasagem (lags)")

    fig.suptitle("Assinaturas típicas: como ler ACF e PACF para identificar ordens p e q",
                 fontsize=12.5, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{OUT}/fig2_acf_pacf.png", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# FIGURA 3 — Validação: origem rolante + erro por horizonte (ilustrativo)
# ----------------------------------------------------------------------
def figura_validacao():
    fig = plt.figure(figsize=(9.8, 6))

    # (a) Esquema de origem rolante
    ax1 = fig.add_axes([0.08, 0.56, 0.88, 0.34])
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 10)
    ax1.axis("off")
    ax1.set_title("(a) Validação por origem rolante (rolling origin / backtesting)",
                  fontsize=11.5, fontweight="bold", loc="left")

    segmentos = [(0, 18, 2), (10, 18, 3), (20, 18, 4), (30, 18, 5), (40, 18, 6)]
    for i, (inicio, largura, y0) in enumerate(segmentos):
        ax1.add_patch(Rectangle((inicio, y0), largura, 0.9, color=AZUL, alpha=0.85))
        ax1.add_patch(Rectangle((inicio + largura, y0), 3, 0.9, color="#c0392b"))
        ax1.text(inicio + largura + 4, y0 + 0.45, "teste", va="center", fontsize=8, color="#c0392b")
    ax1.text(1, 0.6, "treino cresce →\ntestagens sempre no futuro", fontsize=9, color="#333333")
    ax1.annotate("", xy=(78, 8.2), xytext=(64, 8.2),
                 arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))
    ax1.text(64, 8.9, "o modelo é refitado a cada origem", fontsize=9, color="#333333")
    ax1.text(80, 7.9, "erro calculado\napenas no ponto vermelho", fontsize=8.5, color="#333333")

    # (b) Erro médio por horizonte — ilustração didática
    ax2 = fig.add_axes([0.08, 0.08, 0.88, 0.38])
    h = np.arange(1, 13)
    naive = 8 + 2.4 * np.sqrt(h)          # passeio aleatório: erro cresce com √h
    model = 5.2 + 0.9 * np.sqrt(h) + 0.45 * h  # modelo com viés leve
    ax2.plot(h, naive, color=LARANJA, lw=2, marker="o", ms=4, label="Passeio aleatório (naïve)")
    ax2.plot(h, model, color=AZUL, lw=2, marker="s", ms=4, label="Modelo ajustado")
    ax2.set_xlabel("horizonte de previsão (passos à frente)")
    ax2.set_ylabel("erro absoluto médio (MAE) — valores ilustrativos")
    ax2.set_title("(b) Como o erro muda com o horizonte — escolher h pela decisão, não pelo menor erro global",
                  fontsize=11.5, fontweight="bold", loc="left")
    ax2.legend(frameon=False, loc="upper left")
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{int(v)}"))

    fig.savefig(f"{OUT}/fig3_validacao.png", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# FIGURA 4 — Dados reais do BCB: IPCA (mensal e 12 meses) e Meta Selic
# ----------------------------------------------------------------------
IPCA_MESES = ["2024-10", "2024-11", "2024-12",
              "2025-01", "2025-02", "2025-03", "2025-04", "2025-05", "2025-06",
              "2025-07", "2025-08", "2025-09", "2025-10", "2025-11", "2025-12",
              "2026-01", "2026-02", "2026-03", "2026-04", "2026-05", "2026-06",
              "2026-07", "2026-08", "2026-09"]
IPCA_VAL = [0.56, 0.39, 0.52, 0.16, 1.31, 0.56, 0.43, 0.26, 0.24, 0.26, -0.11,
            0.48, 0.09, 0.18, 0.33, 0.33, 0.70, 0.88, 0.67, 0.58, 0.16, 0.07, -0.32, 0.82]

SELIC_PASSOS = [  # (data ISO, meta % a.a.) — mudanças definidas pelo Copom
    ("2023-01-01", 13.75), ("2023-08-03", 13.25), ("2023-09-21", 12.75),
    ("2023-11-02", 12.25), ("2023-12-14", 11.75), ("2024-02-01", 11.25),
    ("2024-03-21", 10.75), ("2024-05-09", 10.50), ("2024-09-19", 10.75),
    ("2024-11-07", 11.25), ("2024-12-12", 12.25), ("2025-01-30", 13.25),
    ("2025-03-20", 14.25), ("2025-05-08", 14.75), ("2025-06-19", 15.00),
    ("2026-03-19", 14.75), ("2026-04-30", 14.50), ("2026-06-18", 14.25),
    ("2026-08-06", 14.00), ("2026-09-17", 13.75),
]


def figura_brasil():
    mes_idx = np.arange(len(IPCA_MESES))
    y_mom = np.array(IPCA_VAL)

    # Acumulado 12 meses por composição multiplicativa
    y12 = []
    for i in range(len(y_mom)):
        if i >= 11:
            comp = np.prod(1 + y_mom[i - 11 : i + 1] / 100) - 1
            y12.append(comp * 100)
        else:
            y12.append(np.nan)
    y12 = np.array(y12)

    fig = plt.figure(figsize=(10, 7.6))
    ax1 = fig.add_axes([0.07, 0.53, 0.89, 0.40])
    cores = ["#c0392b" if v < 0 else AZUL for v in y_mom]
    ax1.bar(mes_idx, y_mom, color=cores, width=0.62)
    ax1.set_xticks(mes_idx)
    labels = [m.replace("-", "/") for m in IPCA_MESES]
    ax1.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax1.set_ylabel("variação mensal (%)")
    ax1.set_title("(a) IPCA — variação mensal (%) | fonte: BCB, série SGS 433 (dados reais)",
                  fontsize=11.5, fontweight="bold", loc="left")
    ax1.axhline(0, color="#333333", lw=0.8)
    ax1.yaxis.set_major_formatter(FuncFormatter(fmt_br))

    ax1b = ax1.twinx()
    ax1b.plot(mes_idx, y12, color=LARANJA, lw=2.2, marker="D", ms=4,
              label="acumulado em 12 meses (%)")
    ax1b.set_ylabel("acumulado 12 meses (%)", color=LARANJA)
    ax1b.tick_params(axis="y", colors=LARANJA)
    ax1b.yaxis.set_major_formatter(FuncFormatter(fmt_br))
    ax1b.legend(frameon=False, loc="lower left")
    ax1b.set_ylim(0, max(np.nanmax(y12) * 1.25, 8))

    # Painel inferior: Meta Selic (escada)
    ax2 = fig.add_axes([0.07, 0.07, 0.89, 0.34])
    datas = [np.datetime_data(np.datetime64(d[0])) for d in SELIC_PASSOS]
    vals = [d[1] for d in SELIC_PASSOS]
    dt = np.array([np.datetime64(d[0]) for d in SELIC_PASSOS])
    ax2.step(dt, vals, where="post", color=VERDE, lw=2.2)
    ax2.scatter(dt, vals, color=VERDE, s=18, zorder=3)
    ax2.set_ylabel("meta Selic (% ao ano)")
    ax2.set_title("(b) Meta Selic definida pelo Copom (% a.a., escada) | fonte: BCB, série SGS 432 (dados reais)",
                  fontsize=11.5, fontweight="bold", loc="left")
    ax2.grid(axis="y", alpha=0.25)
    import datetime as _dt
    ticks = [_dt.date(2023, 1, 1), _dt.date(2023, 7, 1), _dt.date(2024, 1, 1),
             _dt.date(2024, 7, 1), _dt.date(2025, 1, 1), _dt.date(2025, 7, 1),
             _dt.date(2026, 1, 1), _dt.date(2026, 7, 1)]
    ax2.set_xticks([np.datetime64(t) for t in ticks])
    ax2.set_xticklabels(["jan/23", "jul/23", "jan/24", "jul/24", "jan/25", "jul/25",
                         "jan/26", "jul/26"], fontsize=8.5)
    for x0, y0 in zip(dt, vals):
        if y0 == max(vals) or y0 == 13.75:
            ax2.annotate(f"{y0:.2f}".replace(".", ","), (x0, y0),
                         textcoords="offset points", xytext=(0, 7), fontsize=8,
                         ha="center", color="#333333")

    fig.savefig(f"{OUT}/fig4_brasil.png", bbox_inches="tight")
    plt.close(fig)

    return {
        "ipca_set_2026_mom": y_mom[-1],
        "ipca_12m_set_2026": float(y12[-1]),
        "ipca_12m_computed_check": float(np.prod(1 + y_mom[-12:] / 100) - 1) * 100,
        "selic_atual": SELIC_PASSOS[-1][1],
        "selic_pico": max(v for _, v in SELIC_PASSOS),
    }


if __name__ == "__main__":
    figura_componentes()
    figura_acf_pacf()
    figura_validacao()
    info = figura_brasil()
    with open(f"{OUT}/../numeros_chave.json", "w") as f:
        json.dump(info, f, ensure_ascii=False, indent=2)
    print("Figuras geradas em:", OUT)
    print(json.dumps(info, ensure_ascii=False, indent=2))
