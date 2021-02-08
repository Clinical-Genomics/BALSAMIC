# Provinformation

ProvID
: {{ meta.case_name }}

Analyspanel
: {{ meta.gene_panel_name }} 

Analyskod
: {{ meta.apptag }} 

Analysdatum
: {{ meta.config_date }}

Bioinformatisk analys
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

**Mediansekvensdjup**: Median av sekvenseringsdjup över alla baser inkluderade i analyspanelen.

**Fold 80 base penalty**: Jämnhet av täckningsgraden över alla gener i analyspanlen. Ett värde mellan 1.0-1.8 visar god jämnhet.

**Fragmentlängd, medel:** Medelstorlek av provbiblioteken som laddats på sekvenseringsinstrument. <200bp kan tyda på degraderade provmaterial (t.ex. FFPE), innan biblioteksberedning.  

**Täckningsgrad:** Andel baser som är sekvenserade med ett djup över en specificerad gräns, t.ex. 100X.
