'''
Contains constants and models for analysis or filtering
'''

VARCALL = {
    'tumor_only': {
        'vardict': {
            'filter': 'AF>0',
            'post_query_script': None
        },
        'mutect2': {
            'filter': 'AF>0',
            'post_query_script': None
        }
    },
    'tumor_normal': {
        'strelka': {
            'filter': 'AF>0',
            'post_query_script': None
        },
        'vardict': {
            'filter': 'AF>0',
            'post_query_script': None
        },
        'mutect2': {
            'filter': 'AF>0',
            'post_query_script': None
        }
    }
}
