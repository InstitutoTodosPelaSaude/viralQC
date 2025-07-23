## Install

### Dependencies

```bash
micromamba env create -f env.yml
micromamba activate viralQA
```

### viralQA

```bash
pip install .
```

### Check installation (CLI)

```bash
vqa --help
```

## Usage (CLI)

### get-nextclade-datasets

This command configures local datasets using nextclade.

```bash
vqa get-nextclade-datasets --cores 2
```

## Usage (API)

```bash
vqa-server
```

Go to `http://127.0.0.1:8000/docs`

### Development

Install development dependencies and run `black` into `viral` directory.

```bash
pip install -e ".[dev]"
black viralqa
```
