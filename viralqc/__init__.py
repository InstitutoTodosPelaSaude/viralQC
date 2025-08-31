from importlib.resources import files

PKG_PATH = files("viralqc")
DATASETS_CONFIG_PATH = PKG_PATH.joinpath("config/nextclade_datasets.yml")
GET_NC_PUBLIC_DATASETS_SNK_PATH = PKG_PATH.joinpath("scripts/get_public_datasets.smk")
GET_BLAST_DB_SNK_PATH = PKG_PATH.joinpath("scripts/get_blast_database.smk")
RUN_NEXTCLADE_SNK_PATH = PKG_PATH.joinpath("scripts/run_nextclade.smk")
