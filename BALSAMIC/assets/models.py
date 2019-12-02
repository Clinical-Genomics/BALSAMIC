'''
Contains constants and models for analysis or filtering
'''

# A dictionary of various filters for variant callers
# 'filter': the most basic filters of them all
# 'somatic_filter': a more sophisticated filter to filter out possible non somatic filters
VARCALL = {
    'tumor_only': {
        'vardict': {
            'filter': 'AF>0',
            'somatic_filter': '',
            'post_query_script': None
        },
        'mutect2': {
            'filter': 'AF>0',
            'somatic_filter': '',
            'post_query_script': None
        }
    },
    'tumor_normal': {
        'strelka': {
            'filter': 'AF>0',
            'somatic_filter': '',
            'post_query_script': None
        },
        'vardict': {
            'filter': 'AF>0',
            'somatic_filter': '',
            'post_query_script': None
        },
        'mutect2': {
            'filter': 'AF>0',
            'somatic_filter': '',
            'post_query_script': None
        }
    }
}
