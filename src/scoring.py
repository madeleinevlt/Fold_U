def ss_score(ali):
    score = 0
    L = 0
    for ind, res_q in enumerate(ali.query.residues):
        res_t = ali.template.residues[ind]
        if res_q.ss == "-" or res_t.ss == "-":
            continue
        if res_q.ss_confidence > 7:
            L += 1
            if res_q.ss == res_t.ss:
                pass
            else:
                # print(res_q.ss, res_t.ss)
                score += 1
    score = round((score/L),2) # Revoir la normalisation
    # print(ali.num, score)
    return score


def blossum_score(ali, BLOSSUM_DICT):
    query = ali.query.residues
    template = ali.template.residues
    score_ali = 0
    size = len(query)
    for i in range(size):
        if query[i].name == "-" or template[i].name == "-":
            continue
        else:
            score_ali += BLOSSUM_DICT[query[i].name+template[i].name]
    return(score_ali)
