import os
import zipfile
import glob

import pandas as pd
from Bio import SeqIO

fasta = "ZMays/fasta"
tes = "ZMays/tes"
groups = "ZMays/Agrup"

# Iterate through the zip files
for i in range(1, 11):
	zip_filename = f"{fasta}/chromosome{i}.zip"

	#Check if the zip file exists
	if os.path.exists(zip_filename):
		try:
			with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
				zip_ref.extractall()
				print(f"Extração bem-sucedida {zip_filename}")
		except zipfile.BadZipFile:
			print(f"Erro: {zip_filename} não é um arquivo zip.")
		except Exception as e:
			print(f"Ocorreu um erro durante a extração {zip_filename}: {e}")
	else:
		print(f"Warning: {zip_filename} não encontrado.")

# Use glob to find all .zip files
zip_files = glob.glob("*.zip")

# Iterate through the found zip files
for zip_filename in zip_files:
	try:
		os.remove(zip_filename)
		print(f"Removido com sucesso {zip_filename}")
	except OSError as e:
		print(f"Erro: {zip_filename}: {e}")
		
file_path = (tes + "/TEAnnotationFinal_Helitron.gff3")

with open(file_path, "r") as file:
	# Use list comprehension para pegar as primeiras 10 linhas
	for i, line in enumerate(file):
		print(line.strip())  # strip() remove espaços extras e quebras de linha duplicadas
		if i == 9:  # Exibe apenas as primeiras 10 linhas
			break
		
def add_gff3_to_dataframe(df, file_path):
	"""
	Lê um arquivo GFF3, extrai os dados e adiciona ao DataFrame existente.

	Parâmetros:
		df (pd.DataFrame): DataFrame existente para o qual os dados serão adicionados.
		file_path (str): Caminho do arquivo GFF3.

	Retorna:
		pd.DataFrame: DataFrame atualizado com os novos dados.
	"""
	# Definir os nomes das colunas do arquivo GFF3
	column_names = ["Chr", "SourceAnnotation", "COS", "Start", "End", "Score", "Strand", "Phase", "Attributes"]

	# Lista para armazenar os dados do novo arquivo
	data = []

	try:
		with open(file_path, "r") as file:
			for line in file:
				if not line.startswith("#"):  # Ignorar linhas de comentários
					parts = line.strip().split("\t")
					if len(parts) == len(column_names):  # Verificar o número correto de colunas
						data.append(parts)
					else:
						print(f"Warning: Skipping line with incorrect number of columns: {line.strip()}")

		# Criar um DataFrame com os dados lidos
		new_df = pd.DataFrame(data, columns=column_names)

		# Adicionar os novos dados ao DataFrame existente
		updated_df = pd.concat([df, new_df], ignore_index=True)
		return updated_df

	except FileNotFoundError:
		print(f"Error: File not found at {file_path}")
	except Exception as e:
		print(f"An error occurred: {e}")

	return df  # Retorna o DataFrame original caso haja um erro

helitron = "ZMays/tes/TEAnnotationFinal_Helitron.gff3"
line = "ZMays/tes/TEAnnotationFinal_LINE.gff3"
ltr = "ZMays/tes/TEAnnotationFinal_LTR.gff3"
mite = "ZMays/tes/TEAnnotationFinal_MITE.gff3"
sine = "ZMays/tes/TEAnnotationFinal_SINE.gff3"
tir = "ZMays/tes/TEAnnotationFinal_TIR.gff3"

df = ""
df = pd.DataFrame(columns=["Chr", "SourceAnnotation", "COS", "Start", "End", "Score", "Strand", "Phase", "Attributes"])

df = add_gff3_to_dataframe(df, helitron)
df = add_gff3_to_dataframe(df, line)
df = add_gff3_to_dataframe(df, ltr)
df = add_gff3_to_dataframe(df, mite)
df = add_gff3_to_dataframe(df, sine)
df = add_gff3_to_dataframe(df, tir)

print(df.head())

#df = df[df['Strand'] == '+']
df.reset_index(drop=True, inplace=True)

print(df.head())
print(f"\n\nNúmero de entradas: {len(df)}")

# Ordernar o DataFrame
df['Start'] = df['Start'].astype(int)
df['End'] = df['End'].astype(int)
df_sorted = df.sort_values(by=['Chr', 'Start'])

# Reset the index of the sorted DataFrame
df_sorted.reset_index(drop=True, inplace=True)

print(df_sorted.head())
print(f"\n\nNúmero de entradas: {len(df_sorted)}")

# Lista dos tipos permitidos
types = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]

# Filtrar o DataFrame para manter apenas os tipos permitidos
df = df[df["Chr"].isin(types)]
df.reset_index(drop=True, inplace=True)

# Sumarização
summary = df["Chr"].value_counts()

for index, value in summary.items():
	print(f"{index}: {value}")

print(f"\n\nNúmero de entradas: {len(df)}")

dfs_metrics = []

i = 0
for f in glob.glob(groups + "/*.csv"):
	df_metric = pd.read_csv(f)
	df_metric = df_metric.sort_values(by='nameseq')

	columns = [f"{f.replace('ZeaMays-TEsClassification/ZMays/Agrup/', '').replace('.csv', '')}_{column}" for column in df_metric.columns]
	if i != 0:
		df_metric.drop(columns=['nameseq'], inplace=True)
		columns.remove(columns[0])
	else:
		columns[0] = 'nameseq'

	df_metric.columns = columns

	df_metric.reset_index(drop=True, inplace=True)
	dfs_metrics.append(df_metric)
	i += 1

df_metrics = pd.concat(dfs_metrics, axis=1)

mask = df_metrics['nameseq'].apply(lambda x: isinstance(x, str))
df_filtrado = df_metrics.loc[mask]

mascara = df_metrics['nameseq'].map(type) == str

# Aplica a máscara para filtrar o DataFrame
df_metrics_filtrado = df_metrics[mascara]

print("DataFrame original:")
print(df_metrics)
print("\nDataFrame filtrado (apenas strings em 'nameseq'):")
print(df_metrics_filtrado)

def juntar_csvs_em_dataframe(arquivo_zip1: str, arquivo_zip2: str, arquivo_zip3: str, arquivo_zip4: str, arquivo_zip5: str) -> pd.DataFrame | None:
	try:
		#with zipfile.ZipFile(arquivo_zip1, 'r') as zipfile1:
		#  with zipfile.ZipFile(arquivo_zip2, 'r') as zipfile2:
		#    with zipfile.ZipFile(arquivo_zip3, 'r') as zipfile3:
		#      with zipfile.ZipFile(arquivo_zip4, 'r') as zipfile4:
		#        with zipfile.ZipFile(arquivo_zip5, 'r') as zipfile5:
		#zipfile1.extractall()
		#rint(f"Extração bem-sucedida {arquivo_zip1}")
		#ipfile2.extractall()
		#print(f"Extração bem-sucedida {arquivo_zip2}")
		#zipfile3.extractall()
		#print(f"Extração bem-sucedida {arquivo_zip3}")
		#print(f"Extração bem-sucedida {arquivo_zip4}")
		#zipfile5.extractall()
		#print(f"Extração bem-sucedida {arquivo_zip5}")

		arquivo_csv1 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_1_de_5.csv'
		arquivo_csv2 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_2_de_5.csv'
		arquivo_csv3 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_3_de_5.csv'
		arquivo_csv4 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_4_de_5.csv'
		arquivo_csv5 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_5_de_5.csv'

		try:
			print(f"\nLendo o primeiro arquivo para junção: '{arquivo_csv1}'...")
			df1 = pd.read_csv(arquivo_csv1)
			# Seleciona apenas as colunas cujo nome tem tamanho 6 e começa com A, T, C ou G
			colunas_kmer = [col for col in df1.columns if all(c in {'A', 'T', 'C', 'G'} for c in col[-6:])]
			df1.drop(columns=colunas_kmer, inplace=True)
			print(f"Dimensões de '{os.path.basename(arquivo_csv1)}': {df1.shape[0]} linhas, {df1.shape[1]} colunas")

			print(f"Lendo o segundo arquivo para junção: '{arquivo_csv2}'...")
			df2 = pd.read_csv(arquivo_csv2)
			# Seleciona apenas as colunas cujo nome tem tamanho 6 e começa com A, T, C ou G
			colunas_kmer = [col for col in df2.columns if all(c in {'A', 'T', 'C', 'G'} for c in col[-6:])]
			df2.drop(columns=colunas_kmer, inplace=True)
			print(f"Dimensões de '{os.path.basename(arquivo_csv2)}': {df2.shape[0]} linhas, {df2.shape[1]} colunas")

			print(f"\nLendo o primeiro arquivo para junção: '{arquivo_csv3}'...")
			df3 = pd.read_csv(arquivo_csv3)
			# Seleciona apenas as colunas cujo nome tem tamanho 6 e começa com A, T, C ou G
			colunas_kmer = [col for col in df3.columns if all(c in {'A', 'T', 'C', 'G'} for c in col[-6:])]
			df3.drop(columns=colunas_kmer, inplace=True)
			print(f"Dimensões de '{os.path.basename(arquivo_csv3)}': {df3.shape[0]} linhas, {df3.shape[1]} colunas")

			print(f"Lendo o segundo arquivo para junção: '{arquivo_csv4}'...")
			df4 = pd.read_csv(arquivo_csv4)
			# Seleciona apenas as colunas cujo nome tem tamanho 6 e começa com A, T, C ou G
			colunas_kmer = [col for col in df4.columns if all(c in {'A', 'T', 'C', 'G'} for c in col[-6:])]
			df4.drop(columns=colunas_kmer, inplace=True)
			print(f"Dimensões de '{os.path.basename(arquivo_csv4)}': {df4.shape[0]} linhas, {df4.shape[1]} colunas")

			print(f"\nLendo o primeiro arquivo para junção: '{arquivo_csv5}'...")
			df5 = pd.read_csv(arquivo_csv5)
			# Seleciona apenas as colunas cujo nome tem tamanho 6 e começa com A, T, C ou G
			colunas_kmer = [col for col in df5.columns if all(c in {'A', 'T', 'C', 'G'} for c in col[-6:])]
			df5.drop(columns=colunas_kmer, inplace=True)
			print(f"Dimensões de '{os.path.basename(arquivo_csv5)}': {df5.shape[0]} linhas, {df5.shape[1]} colunas")

			df_concatenado = pd.concat([df1, df2, df3, df4, df5], ignore_index=True)
			print(f"Dimensões do DataFrame final concatenado: {df_concatenado.shape[0]} linhas, {df_concatenado.shape[1]} colunas")
			print(f"Arquivos '{arquivo_csv1}' e '{arquivo_csv2}' juntados em um DataFrame com sucesso!")
			return df_concatenado
		except FileNotFoundError:
			print(f"Erro: Um ou todos os arquivos CSV ('{arquivo_csv1}', '{arquivo_csv2}', '{arquivo_csv3}', '{arquivo_csv4}', '{arquivo_csv5}') não foram encontrados.")
		return None

	#except zipfile.BadZipFile:
	  #print(f"Erro: não é um arquivo zip.")
	except Exception as e:
	  print(f"Ocorreu um erro durante a extração: {e}")


def juntar_dataframes_inner(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame | None:

	if 'nameseq' not in df1.columns or 'nameseq' not in df2.columns:
		print("Erro: A coluna 'nameseq' não está presente em ambos os DataFrames.")
		return None

	df_juntado = pd.merge(df1, df2, on='nameseq', how='inner')
	print(f"DataFrames juntados com 'inner join' na coluna 'nameseq'. Dimensões do resultado: {df_juntado.shape}")
	return df_juntado

arquivo_csv_parte1 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_1_de_5.zip'
arquivo_csv_parte2 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_2_de_5.zip'
arquivo_csv_parte3 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_3_de_5.zip'
arquivo_csv_parte4 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_4_de_5.zip'
arquivo_csv_parte5 = 'ZMays/Agrup/k-mer/divididos/k-mer_agrupado_parte_5_de_5.zip'



df_kmer = juntar_csvs_em_dataframe(arquivo_csv_parte1, arquivo_csv_parte2, arquivo_csv_parte3, arquivo_csv_parte4, arquivo_csv_parte5)

# Use glob to find all .zip files
zip_files = glob.glob("*.zip")

# Iterate through the found zip files
for zip_filename in zip_files:
	try:
		os.remove(zip_filename)
		print(f"Removido com sucesso {zip_filename}")
	except OSError as e:
		print(f"Erro: {zip_filename}: {e}")

df_kmer.head()

df_all_metrics = juntar_dataframes_inner(df_metrics_filtrado, df_kmer)

#df_all_metrics = df_metrics_filtrado.copy()

print(df_all_metrics.head())

print(df_all_metrics.shape)

def join_metrics_inner(df_main, df_metrics_values):

	df_metrics = df_metrics_values.copy()

	# Dividir 'nameseq' em colunas chave
	try:

		split_data = df_metrics['nameseq'].str.split('_', expand=True, n=2)

		if split_data.shape[1] < 3:
			raise ValueError(
				"A coluna 'nameseq' em df_metrics_values não pôde ser dividida em 3 partes (Chr, Start, End) consistentemente."
			)

		df_metrics['_Chr_key_metric'] = split_data[0].astype(str)
		df_metrics['_Start_key_metric'] = split_data[1].astype(str)
		df_metrics['_End_key_metric'] = split_data[2].astype(str)

	except Exception as e:

		print(f"Erro ao processar a coluna 'nameseq': {e}")
		print("Verifique o formato dos dados em 'nameseq'. Esperado: 'Chr_Start_End'.")

		raise

	# Identificar as colunas de métricas que serão transferidas/atualizadas
	metric_value_cols = [
		col for col in df_metrics.columns
		if col not in ['nameseq', '_Chr_key_metric', '_Start_key_metric', '_End_key_metric']
	]

	# Selecionar apenas as colunas chave e de métricas para a junção
	df_metrics_to_join = df_metrics[
		['_Chr_key_metric', '_Start_key_metric', '_End_key_metric'] + metric_value_cols
	]

	df_left = df_main.copy()

	# Adicionar colunas chave temporárias em df_left, convertendo para string para correspondência
	df_left['_Chr_key_df'] = df_left['Chr'].astype(str)
	df_left['_Start_key_df'] = df_left['Start'].astype(str)
	df_left['_End_key_df'] = df_left['End'].astype(str)

	merged_df = pd.merge(
		df_left,
		df_metrics_to_join,
		left_on=['_Chr_key_df', '_Start_key_df', '_End_key_df'],
		right_on=['_Chr_key_metric', '_Start_key_metric', '_End_key_metric'],
		how='inner', # Especifica a junção interna
		suffixes=('', '_newmetric')
	)


	for metric_col_name in metric_value_cols:

		suffixed_metric_col_name = metric_col_name + '_newmetric'

		if suffixed_metric_col_name in merged_df.columns:

			merged_df[metric_col_name] = merged_df[suffixed_metric_col_name]
			merged_df.drop(columns=[suffixed_metric_col_name], inplace=True)

	cols_to_drop = [
		'_Chr_key_df', '_Start_key_df', '_End_key_df',
		'_Chr_key_metric', '_Start_key_metric', '_End_key_metric' # Remover se ainda existirem
	]

	existing_cols_to_drop = [col for col in cols_to_drop if col in merged_df.columns]
	merged_df.drop(columns=existing_cols_to_drop, inplace=True, errors='ignore')

	return merged_df

df_resultado_juncao = join_metrics_inner(df, df_all_metrics)
print(df_resultado_juncao)

df_met_seq = 'df_metricas_seq.csv'

try:
	df_resultado_juncao.to_csv(df_met_seq, index=False, encoding='utf-8')
	print(f"DataFrame salvo com sucesso como '{df_met_seq}'!")
except Exception as e:
	print(f"Ocorreu um erro ao salvar o DataFrame: {e}")