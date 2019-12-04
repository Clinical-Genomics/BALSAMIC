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
            'filter': 'DP',
            'exclude_filter': ('INFO/DP<50 && INFO/AF<0.005 && ((MQ < 55.0 && NM > 1.0) ',
                               '|| (MQ < 60.0 && NM > 2.0)) && INFO/VD<10'),
            'post_query_script': None
        },
        'mutect2': {
            'filter': 'AF>0',
            'exclude_filter': '',
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
