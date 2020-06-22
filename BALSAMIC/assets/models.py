'''
Contains constants and models for analysis or filtering
'''

# A dictionary of various filters for variant callers
# 'filter': the most basic filters of them all
# 'exclude_filter': exclude criteria for variants
# 'somatic_filter': a more sophisticated filter to filter out possible non somatic filters
NGS_FILTER = {
    'tumor_only': {
        'vardict': {
            'general': {
                'AD': ['20', 'balsamic_low_tumor_ad'],
                'MQ': ['50', 'balsamic_low_mq'],
                'DP': ['100', 'balsamic_low_tumor_dp'],
                'AF_min': ['0.05', 'balsamic_low_af'],
                'AF_max': ['1', 'balsamic_af_equal_one']
            }
        }
    },
    'tumor_normal': {
        'vardict': {
            'general': {
                'AD': ['20', 'balsamic_low_tumor_ad'],
                'MQ': ['50', 'balsamic_low_mq'],
                'DP': ['100', 'balsamic_low_tumor_dp'],
                'AF_min': ['0.05', 'balsamic_low_af'],
                'AF_max': ['1', 'balsamic_af_equal_one']
                'TN_af_ration': ['5', 'balsamic_low_tn_af_ratio']
            }
        }
    }
}
