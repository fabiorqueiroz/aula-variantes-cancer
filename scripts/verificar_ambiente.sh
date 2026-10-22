#!/usr/bin/env bash
# Confere se todas as ferramentas da aula estão disponíveis no terminal.
VERDE='\033[0;32m'; VERMELHO='\033[1;31m'; SEM='\033[0m'
falhou=0
checar() {
    local nome="$1"; shift
    if command -v "$nome" >/dev/null 2>&1; then
        echo -e "  ${VERDE}✔${SEM} $(printf '%-12s' "$nome") $("$@" 2>&1 | grep -v JAVA_TOOL_OPTIONS | head -1)"
    else
        echo -e "  ${VERMELHO}✘ $nome não encontrado${SEM}"; falhou=1
    fi
}
echo "Ferramentas da aula:"
checar fastqc      fastqc --version
checar trimmomatic trimmomatic -version
checar bwa         bash -c "bwa 2>&1 | grep Version"
checar samtools    bash -c "samtools --version | head -1"
checar bcftools    bash -c "bcftools --version | head -1"
checar picard      bash -c "picard MarkDuplicates --version 2>&1 | tail -1"
checar Rscript     Rscript -e 'suppressWarnings(suppressMessages(library(maftools))); cat("R", as.character(getRversion()), "· maftools", as.character(packageVersion("maftools")))'
echo
if [ $falhou -eq 0 ]; then
    echo -e "${VERDE}Tudo certo! Pode começar a aula (veja o GUIA_DA_AULA.md).${SEM}"
else
    echo -e "${VERMELHO}Algo está faltando.${SEM} Tente:  micromamba activate aula-variantes"
    echo "Se não resolver, rode novamente:  bash instalar.sh"
    exit 1
fi
