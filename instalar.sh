#!/usr/bin/env bash
# =============================================================================
#  Instala o ambiente da aula (FastQC, Trimmomatic, BWA, SAMtools, Picard,
#  BCFtools, R + maftools) usando o micromamba — um conda "de bolso", que não
#  interfere em nenhum conda/Anaconda que já exista na máquina.
#
#  No GitHub Codespaces isto roda sozinho ao criar o ambiente.
#  Em Linux, WSL (Windows) ou macOS:   bash instalar.sh
# =============================================================================
set -euo pipefail

DIR_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MM_VERSAO=2.4.0
MM_BIN="$HOME/.local/bin/micromamba"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ENV=aula-variantes

echo "==> Instalando o ambiente '$ENV' (leva de 3 a 8 minutos)..."

# ---------------------------------------------------------- 1. micromamba
if [ ! -x "$MM_BIN" ]; then
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64)          PLAT=linux-64 ;;
        Linux-aarch64)         PLAT=linux-aarch64 ;;
        Darwin-x86_64)         PLAT=osx-64 ;;
        Darwin-arm64)          PLAT=osx-arm64 ;;
        *) echo "Sistema não suportado: $(uname -s) $(uname -m)"; exit 1 ;;
    esac
    if ! command -v bzip2 >/dev/null 2>&1 && command -v apt-get >/dev/null 2>&1; then
        echo "==> Instalando bzip2..."
        (sudo apt-get update -qq && sudo apt-get install -y -qq bzip2) >/dev/null 2>&1 || true
    fi
    echo "==> Baixando micromamba $MM_VERSAO ($PLAT)..."
    mkdir -p "$HOME/.local/bin"
    TMP=$(mktemp -d)
    curl -fsSL "https://conda.anaconda.org/conda-forge/$PLAT/micromamba-$MM_VERSAO-0.tar.bz2" -o "$TMP/mm.tar.bz2" \
      || curl -fsSL "https://micro.mamba.pm/api/micromamba/$PLAT/latest" -o "$TMP/mm.tar.bz2"
    tar -xjf "$TMP/mm.tar.bz2" -C "$TMP" bin/micromamba
    mv "$TMP/bin/micromamba" "$MM_BIN"
    rm -rf "$TMP"
fi

# ----------------------------------------------- 2. ambiente com as ferramentas
PLAT_ENV=()
if [ "$(uname -s)" = "Darwin" ] && [ "$(uname -m)" = "arm64" ]; then
    # Nem todas as ferramentas do bioconda têm versão nativa para Apple Silicon:
    # usamos a versão Intel (roda via Rosetta 2).
    PLAT_ENV=(--platform osx-64)
fi
if "$MM_BIN" env list 2>/dev/null | grep -q "envs/$ENV\$"; then
    echo "==> Ambiente '$ENV' já existe — atualizando, se necessário."
    "$MM_BIN" install -y -n "$ENV" -f "$DIR_REPO/environment.yml" "${PLAT_ENV[@]}" >/dev/null
else
    "$MM_BIN" create -y -f "$DIR_REPO/environment.yml" "${PLAT_ENV[@]}"
fi
"$MM_BIN" clean -a -y >/dev/null 2>&1 || true

# ------------------------------- 3. ativar o ambiente automaticamente no terminal
touch "$HOME/.bashrc"
for RC in "$HOME/.bashrc" "$HOME/.zshrc"; do
    [ -f "$RC" ] || continue
    if ! grep -q ">>> aula-variantes >>>" "$RC"; then
        SHELL_NOME=$(basename "$RC" | sed 's/^\.//; s/rc$//')
        cat >> "$RC" <<EOF

# >>> aula-variantes >>>  (para remover, apague este bloco)
export MAMBA_ROOT_PREFIX="$MAMBA_ROOT_PREFIX"
eval "\$("$MM_BIN" shell hook -s $SHELL_NOME)"
micromamba activate $ENV
# <<< aula-variantes <<<
EOF
    fi
done

echo
echo "============================================================"
echo "  Ambiente pronto!"
echo "  Abra um NOVO terminal (ou rode:  source ~/.bashrc )"
echo "  e confira com:   bash scripts/verificar_ambiente.sh"
echo "============================================================"
