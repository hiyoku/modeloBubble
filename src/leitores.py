"""
    Arquivo contendo os leitoers dos arquivos DAT
    Cada função é lê um tipo de .DAT
"""


def leitor_frequencia_fof2(arquivo_entrada=None):
    """
        Função para ler um arquivo DAT
        Retorna um vetor
    """
    result = []   # Criando um vetor vazio

    try:
        arquivoDAT = open(arquivo_entrada, 'r')   # Abrindo um Arquivo
        for linha in arquivoDAT:
            linha = linha.split("  ")
            result.append(float(linha[1]))
    except Exception as e:
        print(e)
    finally:
        return result


def leitor_campo_eletrico_zonal(arquivo_entrada=None):
    pass