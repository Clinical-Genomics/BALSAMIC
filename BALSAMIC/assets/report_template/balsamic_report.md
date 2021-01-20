# kvalitetsrapport: klinisk sekvensering av cancer prover


##Sammanfattning

{% set table_qc_header = [] %}
{% for qc_name, qc_name_meta in analysis_results.qc.items() -%}
{{ table_qc_header.append(qc_name_meta[lang]) or ''}}
{% endfor -%}

| Prov-id | Provtyp | Familj | Analysdatum | Analyspanel | {{ table_qc_header|join('|') }} |
| --- | --- | --- | --- | --- | --- | --- | --- |
{%- for sample in meta.sample_map.keys() %} 
| {{ meta.sample_map[sample] }} | {{ meta.sample_type[sample] }} | {{ meta.case_name }} | {% if meta.date %} {{ meta.date }} {% else %} {{ meta.now }} {% endif %} | {{ meta.gene_panel_name }} | {{ to_report["qc"][sample]|join('|') }} | 
{%- endfor %}


<br>

## Täckningsgrad för vald analyspanel
{% set table_coverage_header = [] %}
{% for coverage_name, coverage_name_meta in analysis_results.coverage.items() -%}
{{ table_coverage_header.append(coverage_name_meta[lang]) or ''}}
{% endfor %}

| Prov-id | Provtyp | {{ table_coverage_header|join('|') }} |
| --- | --- | --- | --- | --- | --- | --- |
{%- for sample in meta.sample_map.keys() %} 
| {{ meta.sample_map[sample] }} | {{ meta.sample_type[sample] }} | {{ to_report.coverage[sample]|join('|') }} |
{%- endfor %}

### Förklaringar
Värdena presenterade för Mediansekvensdjup och Täckningsgrad är efter borttagning av duplikata läsningar.

**Mediansekvensdjup**: Median av seknvenseringsdjup över alla baser inkluderade i analyspanelen.

**Fold-80**: Jämnhet av täckningsgraden över alla gener i analyspanlen. Ett värde mellan 1.0-1.8 visar god jämnhet.

**Medelinsättningsstorlek:** Medelstorlek av provbiblioteken som laddats på sekvenseringsinstrument. <200bp tyder på degraderade provmaterial (t.ex. FFPE), innan biblioteksberedning.  

**Täckningsgrad:** Andel baser som är sekvenserade med ett djup över en specificerad gräns, t.ex. 100X.
