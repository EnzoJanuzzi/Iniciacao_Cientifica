# Iniciação Científica FAPESP — Simulações DFT
Bolsista: Enzo Januzzi Xavier

Orientador: Prof. James Moraes de Almeida (ILUM/CNPEM)

Período de vigência: 01/11/2025 a 31/10/2026

Este repositório reúne os arquivos relevantes, dados e scripts utilizados durante a Iniciação Científica financiada pela FAPESP, do processo 2025/20105-9 vinculado ao processo 2023/09820-2. O projeto envolve o uso de métodos de Teoria do Funcional da Densidade (DFT) para o estudo de propriedades eletrônicas e estruturais de nanopartículas bimetálicas (Cu@Ru)com potencial aplicação catalítica. A maioria desses códigos gera dezenas de arquivos por pasta, mas apenas os arquivos de saída utilizados foram adicionados para facilitar a reprodutibilidade do projeto. Os outros arquivos podem ser encontrados no repositório nacional Materials Informatics, através da seguinte [coleção](https://www.midb.cloud/shared_collection?key=c8427432f4f44ec0baead49325d9db8a).

## Estrutura do Repositório
- Scripts: Scripts gerais para automatização dos processos e organização dos dados.
- VASP: Arquivos de entrada e saída das simulações realizadas com o software VASP, para obter as estruturas, trajetórias e energias de cada intermediário adsorvido na nanopartícula. Esses dados foram usados para a avaliação da estabilidade dos sistemas e para o cálculo dos ΔG.
- Python_codes: Arquivos utilizados para o pós-processamento e geração dos gráficos de energia livre, contendo os dados extraídos das simulações e o fluxograma para o cálculo e plotagem.
