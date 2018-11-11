# Formalign.eu
Format and display multiple sequence alignments (nucleotide or protein).

This app is hosted at [formalign.eu](https://formalign.eu)

## Usage
```bash
pip install numpy
pip install requirements.txt
pip install uwsgi
uwsgi --ini formalign_uwsgi.ini
```

## Test
```bash
python manage.py test base/tests
```
