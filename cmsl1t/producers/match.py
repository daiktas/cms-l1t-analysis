import math
# TODO: make it into a proper producer


def get_matched_obj_index(obj, objects, deltaR=0.3):
    obj_eta, obj_phi = obj
    closest_index = 9999
    minDeltaR = deltaR

    for i, o in enumerate(objects):
        o_eta, o_phi = o
        dEta = obj_eta - o_eta
        dPhi = obj_phi - o_phi
        dR = math.sqrt(dEta**2 + dPhi**2)
        if dR < minDeltaR:
            minDeltaR = dR
            closest_index = i
    return closest_index


def get_matched_l1_jet(recoJet, l1Jets, deltaR=0.3):
    recoObj = (recoJet['eta'], recoJet['phi'])
    l1Objects = [(j['eta'], j['phi']) for j in l1Jets]
    index = get_matched_obj_index(recoObj, l1Objects, deltaR)
    if index < len(l1Jets):
        return l1Jets[index]
    else:
        return None
