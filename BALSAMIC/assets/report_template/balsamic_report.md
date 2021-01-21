# kvalitetsrapport: klinisk sekvensering av cancer prover

# Provinformation

Familj
: {{ meta.case_name }}

Analysispanel
: {{ meta.gene_panel_name }} 

Analysisdatum
: {{ meta.now }}

Bioinformatikanlys:
: {{ meta.bioinformatic }}

# Sammanfattning

| Prov-id | Provtyp | {{ meta.qc_table_header|join('|') }} |
| :---: | :---: | :---: | :---: | :---: |
{%- for sample, value in meta.qc_table_content.items() %} 
{{ value|join('|') }} |
{%- endfor %}


# Täckningsgrad för vald analyspanel
| Prov-id | Provtyp | {{ meta.coverage_table_header|join('|') }} |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
{%- for sample, value in meta.coverage_table_content.items() %} 
{{ value|join('|') }} |
{%- endfor %}

# Förklaringar
Värdena presenterade för Mediansekvensdjup och Täckningsgrad är efter borttagning av duplikata läsningar.

**Mediansekvensdjup**: Median av seknvenseringsdjup över alla baser inkluderade i analyspanelen.

**Fold-80**: Jämnhet av täckningsgraden över alla gener i analyspanlen. Ett värde mellan 1.0-1.8 visar god jämnhet.

**Medelinsättningsstorlek:** Medelstorlek av provbiblioteken som laddats på sekvenseringsinstrument. <200bp tyder på degraderade provmaterial (t.ex. FFPE), innan biblioteksberedning.  

**Täckningsgrad:** Andel baser som är sekvenserade med ett djup över en specificerad gräns, t.ex. 100X.

**Obs**: BALSAMIC-analys är inte en ackrediterad metod.
