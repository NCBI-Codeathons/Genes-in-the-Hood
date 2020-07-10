# Install requirements

```bash
python3.7 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
```

# Execute application
```bash
cat testdata.csv | python find_genes_by_location.py --packages-dir /usr/local/data/assemblies/packages
```

