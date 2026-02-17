# MHGuide to CSV

Anwendung, um Angaben aus einer MHGuide JSON Datei in eine CSV-Datei zu exportieren und zur Tumordokumentation zu
verwenden.

## Verwendung

Diese Anwendung erzeugt aus einer MHGuide JSON Datei eine CSV Datei.
Dabei wird der Dateiname beibehalten, die Dateiendung jedoch durch `.csv` ersetzt.

Übernommen werden alle Varianten, die als '(Likely) oncogenic' oder '(Likely) benign' markiert sind.

```
Usage: mhguide-to-csv [OPTIONS] <INPUT_FILE>

Arguments:
  <INPUT_FILE>  Zu lesende JSON-Datei

Options:
      --all-variants  Alle Varianten verwenden, nicht nur '(Likely) oncogenic' oder '(Likely) benign'
  -h, --help          Print help
  -V, --version       Print version
```

Sollen alle Varianten verwendet werden, dann kann dies mit `--all-variants` angegeben werden.

## Enthaltene Liste mit Genen

Es ist eine Liste mit über 43000 Genen
von [https://genenames.org](https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_pub_chrom_map&col=md_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit)
enthalten.

Diese Liste der Gene unterliegt der folgenden Lizenz und ist frei
verfügbar: [Creative Commons Public Domain (CC0) License](https://creativecommons.org/public-domain/cc0/).

Achtung! Es werden nur die aktuell anerkannten Symbolnamen verwendet.
