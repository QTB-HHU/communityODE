import numpy as np

#(PA, F, A, P)
#minimal media objective
objectives0bac_mm = {8: 0.00394660755189, 16: 0.0214594672552, 40: 0.133019126626, 48: 0.325143866598, 64: 0.819741886626, 72: 1.01158000125, 104: 1.12067725568, 112: 1.18736790528, 120: 1.19198783372, 128: 1.23312503242, 152: 1.20902767308, 168: 1.10392536294, 176: 1.09602488201, 184: 0.950595341459, 208: 1.01549725646, 216: 0.96464403229, 232: 0.850799396424, 240: 0.869247193626, 248: 0.704145680136, 264: 0.481277258767, 272: 0.504451680826, 288: 0.394302870761}

#Sampling times:  64.,  120.,  176.
#STARTING {'PSEUDOALTERO': 0.85, 'FLAVO': 0.05, 'ALTEROMO': 0.01, 'PSEUDOMO': 0.005}
objectivesNbac_mm = {64:  {'PSEUDOALTERO': 0.20, 'FLAVO': 0.09, 'ALTEROMO': 0.21, 'PSEUDOMO': 0.18},
                     120: {'PSEUDOALTERO': 0.33, 'FLAVO': 0.03, 'ALTEROMO': 0.55, 'PSEUDOMO': 0.03},
                     176: {'PSEUDOALTERO': 0.37, 'FLAVO': 0.02, 'ALTEROMO': 0.18, 'PSEUDOMO': 0.39}}

#full media objective
objectives0bac_fm = {8: 0.0495192526186, 16: 0.0444984512436, 40: 0.162142228696, 48: 0.604924134947, 64: 0.733476718424, 72: 0.918735231869, 104: 1.03685842214, 112: 1.09875658872, 120: 1.13416673961, 128: 1.10830730003, 152: 0.859050588654, 168: 0.82087016665, 176: 0.843908465158, 184: 0.623807803815, 208: 0.68209676795, 216: 0.623958554326, 232: 0.555624957592, 240: 0.535425524594, 248: 0.478181475877, 264: 0.199115838297, 272: 0.281882401222, 288: 0.303430723727}

objectivesNbac_fm = {64:  {'PSEUDOALTERO': 0.07, 'FLAVO': 0.50, 'ALTEROMO': 0.11, 'PSEUDOMO': 0.01},
                     120: {'PSEUDOALTERO': 0.44, 'FLAVO': 0.46, 'ALTEROMO': 0.06, 'PSEUDOMO': 0.01},
                     176: {'PSEUDOALTERO': 0.57, 'FLAVO': 0.08, 'ALTEROMO': 0.18, 'PSEUDOMO': 0.12}}



def getObjective(pseudoaltero, flavo, altero, pseudo):

    bkeys = []
    if pseudoaltero:
        bkeys.append('PSEUDOALTERO')
    if flavo:
        bkeys.append('FLAVO')
    if altero:
        bkeys.append('ALTEROMO')
    if pseudo:
        bkeys.append('PSEUDOMO')

    if len(bkeys) < 1:
        return objectives0bac_fm, objectives0bac_mm

    obj_fm = {}
    obj_mm = {}
    for tk in objectivesNbac_mm.keys():
        fm = []
        mm = []
        for bk in bkeys:
            fm.append(objectivesNbac_fm[tk][bk])
            mm.append(objectivesNbac_mm[tk][bk])
        obj_fm[tk] = doRatios(fm)
        obj_mm[tk] = doRatios(mm)
    return obj_fm, obj_mm

def doRatios(nums):
    n = np.array(nums)
    tot = float(sum(n))
    p = tuple(np.floor(1000*n/tot)/1000)
    return p

