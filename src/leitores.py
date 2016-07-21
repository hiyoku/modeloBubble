"""
    Arquivo contendo os leitoers dos arquivos DAT
    Cada funcao le um tipo de .DAT
"""
import numpy as np


def leitor_frequencia_fof2(arquivo_entrada=None):
    """
        Leitor de Frequencia foF2 21:45 UT
        Funcao para ler um arquivo DAT
        Retorna um vetor
    """
    result = []   # Criando um vetor vazio

    try:
        arquivoDAT = open(arquivo_entrada, 'r')   # Abrindo um Arquivo

        # Quebrando a linha em vetor de 2 posicoes (XH, F2)
        for linha in arquivoDAT:
            # Quebrando a linha em vetor de 2 posicoes (XH, F2)
            linha = linha.split("  ")

            # Adicionando a posicao 2 em um vetor
            result.append(np.float64(float(linha[1])))

    except Exception as e:
        print("Erro no leitor de frequencia foF2\n" + str(e))

    finally:
        return result


def leitor_campo_eletrico_zonal(arquivo_entrada=None):
    """
        Leitor de Campo Electrico Zonal = BG * dhF / dt
        Funcao para ler um arquivo DAT
        Retorna um vetor
    """
    result = []

    try:
        arquivoDAT = open(arquivo_entrada, 'r')   # Abrindo um Arquivo

        for linha in arquivoDAT:
            # Quebrando a linha em vetor de 2 posicoes (hora, drift)
            linha = linha.split("  ")

            # Adicionando a posicao 2 em um vetor
            result.append(np.float64(float(0.19E-04 * float(linha[1]))))

    except Exception as e:
        print("Erro no leitor de campo electrico zonal\n" + str(e))
    finally:
        return result


def leitor_vento_termosferico(arquivo_entrada=None):
    """
        Leitor de vento termosferico UYZ (zonal), UX (Meridional)
        Funcao para ler um arquivo DAT
        Retorna um vetor
    """
    result = []
    result2 = []

    try:
        arquivoDAT = open(arquivo_entrada, 'r')   # Abrindo um Arquivo

        for linha in arquivoDAT:
            # Quebrando a linha em vetor de 3 posicoes (H1, UM, UZ)
            linha = linha.split("  ")

            result.append(np.float64(float(linha[1])))
            result2.append(np.float64(float(linha[2])))

    except Exception as e:
        print("Erro no leitor de vento termosferico\n" + str(e))
    finally:
        return result, result2
