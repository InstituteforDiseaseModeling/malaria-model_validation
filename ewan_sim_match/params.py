sites = ['6']
exp_name = "model_validation_ewan_eirs3_%s" % sites[0]
nSims = 2
simulation_duration = 365*60

summary_report_age_bin_by_site = {'1': [0.4, 0.9, 2, 3, 4, 9, 14],
                                  '2': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 19, 39, 59, 85],
                                  '3': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 19, 39, 59, 85],
                                  '4': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60, 85],
                                  '5': [4, 9, 14, 29, 44, 85],
                                  '6': [1, 2, 6, 10, 14],
                                  '7': [1, 2, 4, 9, 14, 19, 39, 59, 85],
                                  '8': [1, 2, 4, 9, 14, 19, 39, 59, 85],
                                  '9': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60, 85],
                                  }

summary_report_age_bins = summary_report_age_bin_by_site[sites[0]]




# exp_name = "model_validation_age_inc_sweepSeasonEIRCM"
# nSims = 2
# simulation_duration = 365*60
#
# seasonality_names = ['highSeason', 'midSeason', 'flatSeason']
# eirs = [1, 3, 5, 10, 20, 30, 50, 100, 200, 400]
# cms = [0, 0.35, 0.5]
# sites = []
# for ss in range(len(seasonality_names)):
#     for ee in eirs:
#         for cm in cms:
#             sites = sites + ['%s_%i_%i' % (seasonality_names[ss], ee, cm*100)]
#
# summary_report_age_bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60, 85]



