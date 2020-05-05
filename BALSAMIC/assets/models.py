'''
Contains constants and models for analysis or filtering
'''

# A dictionary of various filters for variant callers
# 'filter': the most basic filters of them all
# 'exclude_filter': exclude criteria for variants
# 'somatic_filter': a more sophisticated filter to filter out possible non somatic filters
VARCALL = {
    'tumor_only': {
        'vardict': {
            'exclude_filter': 'FMT/AD[@tumor_sample_names:1] < 10 || FMT/DP[@tumor_sample_names] < 100',
            'exclude_filter_string': 'balsamic_vardict_softfilter',
            'exclude_filter_mode': '-m "+" -s',
            'post_query_script': None
        },
        'mutect2': {
            'exclude_filter': ('FMT/AD[@tumor_sample_names:1]<10 || sum(FMT/AD)<100 ',
                               '|| FMT/AF[@tumor_sample_names]>0.75',
                               '|| (FMT/QSS[@tumor_sample_names:1]/FMT/AD[@tumor_sample_names:1])<20'),
            'exclude_filter_string': 'balsamic_mutect2_softfilter',
            'exclude_filter_mode': '-m "+" -s',
            'post_query_script': None
        }
    },
    'tumor_normal': {
        'strelka': {
            'filter': 'AF>0',
            'exclude_filter': '',
            'post_query_script': None
        },
        'vardict': {
            'filter': 'AF>0',
            'exclude_filter': '',
            'post_query_script': None
        },
        'mutect2': {
            'filter': 'AF>0',
            'exclude_filter': '',
            'post_query_script': None
        }
    }
}
