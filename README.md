# MHGuide to CSV

Anwendung, um Angaben aus einer MHGuide JSON Datei in eine CSV-Datei zu exportieren und zur Tumordokumentation zu
verwenden.

## Verwendung

Diese Anwendung erzeugt aus einer MHGuide JSON Datei eine CSV Datei.
Dabei wird der Dateiname beibehalten, die Dateiendung jedoch durch `.csv` ersetzt.

Mit dem Parameter `--xlsx` kann die Ausgabe als XLSX-Datei im Format für Excel 2007-365 erfolgen.
Auch hier wird der Dateiname beibehalten, die Dateiendung jedoch durch `.xlsx` ersetzt.

Übernommen werden alle Varianten, sofern nicht anders angegeben, die als '(Likely) oncogenic' markiert oder in
`REPORT_NARRATIVE` aufgeführt sind.

```
Usage: mhguide-to-csv [OPTIONS] <INPUT_FILE>

Arguments:
  <INPUT_FILE>  Zu lesende JSON-Datei

Options:
      --all-variants  Alle Varianten verwenden, nicht nur '(Likely) oncogenic' oder aus 'REPORT_NARRATIVE'
      --oncogenic     Nur Varianten mit '(Likely) oncogenic' verwenden, keine aus 'REPORT_NARRATIVE'
      --no-artifacts  Entferne Artefakte aus 'REPORT_NARRATIVE'
      --xlsx          Exportiere im XLSX-Format (Excel 2007-365)
  -h, --help          Print help
  -V, --version       Print version
```

Sollen alle Varianten verwendet werden, dann kann dies mit `--all-variants` angegeben werden.
Um nur Varianten mit '(Likely) oncogenic' zu verwenden, dann kann dies mit `--oncogenic` angegeben werden.

Mit der Option `--no-artifacts` werden alle Varianten, die als Artefakte gekennzeichnet sind, aus der Liste entfernt.
Dies ist nur möglich, wenn keine der folgenden Optionen verwendet wird: `--all-variants` oder `--oncogenic`.

### RNA Fusionen

RNA Fusionen werden exportiert, wenn die Angaben in der JSON-Datei unter `REPORT_NARRATIVE` vorhanden sind.
Hierbei ist das Format einzuhalten:

* Jede Fusion wird in einer neuen Zeile gelistet.
* Jede Teilangabe der Fusion wird durch ein Semikolon getrennt.

Beispiel:
`ABCD1(ex 1)::ABCD2(ex 2); Transcript ID: NM_012345.4/NM_012456.2; Strand: -; Breakpoint: chr19:12345678/chr19:13456789; Supporting read pairs: 1234`

Daraus werden folgende Angaben exportiert:

| Feld                           | Wert                                      |
|--------------------------------|-------------------------------------------|
| H-Nummer                       | ...                                       | 
| Referenz-Genom                 | ...                                       | 
| Ergebnis                       | RNA Fusion                                | 
| Gen                            | ABCD1                                     | 
| Fusioniertes Gen               | ABCD2                                     | 
| 5' Partner EnsemblID           | ENSG00000101986                           | 
| 5' Partner HGNC ID             | HGNC:61                                   | 
| 5' Partner HGNC Name           | ATP binding cassette subfamily D member 1 | 
| 5' Partner Transcript ID       | NM_012345.4                               | 
| 5' Partner Exon ID             | Exon1                                     | 
| 5' Partner Transcript Position | 12345678                                  | 
| 5' Partner Strand              | -                                         | 
| 3' Partner EnsemblID           | ENSG00000173208                           | 
| 3' Partner HGNC ID             | HGNC:66                                   | 
| 3' Partner HGNC Name           | ATP binding cassette subfamily D member 2 | 
| 3' Partner Transcript ID       | NM_012456.2                               | 
| 3' Partner Exon ID             | Exon2                                     | 
| 3' Partner Transcript Position | 13456789                                  | 
| 3' Partner Strand              | -                                         | 
| Number reported reads          | 1234                                      | 
| Pathogenitätsklasse            |                                           | 

## Enthaltene Liste mit Genen

Es ist eine Liste mit rund 45000 Genen
von [https://genenames.org](https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=gd_app_name&col=gd_pub_chrom_map&col=md_ensembl_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit)
enthalten.

Diese Liste der Gene unterliegt der folgenden Lizenz und ist frei
verfügbar: [Creative Commons Public Domain (CC0) License](https://creativecommons.org/public-domain/cc0/).

Achtung! Es werden nur die aktuell anerkannten Symbolnamen verwendet.
