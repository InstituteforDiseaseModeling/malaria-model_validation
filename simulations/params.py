# note: this setup where we manually specify one site at a time is only temporary and will be replaced with a more
#    automatic commissioning process

# site-matched simulations
# sites = ['chonyi_1999']
# sites = ['ngerenya_1999']
# sites = ['dielmo_1990']
# sites = ['ndiop_1993']
# sites = ['ebolakounou_1997']
# sites = ['koundou_1997']
sites = ['matsari_1970']
# sites = ['rafin_marke_1970']
# sites = ['sugungum_1970']
# sites = ['navrongo_2000']
# sites = ['laye_2007']
# sites = ['dapelogo_2007']
# sites = ['dongubougou_1999']
# sites = ['sotuba_1999']
exp_name = "validation_%s" % sites[0]
nSims = 5


# # sweep over site characteristics
# exp_name = 'test_validation_site_char_sweep'
# nSims = 2
# seasonality_names = ['highSeason', 'midSeason', 'flatSeason']
# eirs = [1, 3, 5, 10, 20, 30, 50, 100, 200, 400]
# cms = [0, 0.35, 0.5]
# sites = []
# for ss in range(len(seasonality_names)):
#     for ee in eirs:
#         for cm in cms:
#             sites = sites + ['%s_%i_%i' % (seasonality_names[ss], ee, cm*100)]
#
# sites = sites[:5]



