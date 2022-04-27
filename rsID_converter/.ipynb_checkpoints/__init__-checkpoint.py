import tabix

GCF_path = '/commons/Reference/dbSNP/b155/GCF_000001405.39.gz'
GCF_H_path = '/commons/Reference/dbSNP/b155/GCF_000001405.39H_sorted.gz'

Chr2RefSeq_accession = {'1':'NC_000001.11','2':'NC_000002.12','3':'NC_000003.12','4':'NC_000004.12','5':'NC_000005.10','6':'NC_000006.12','7':'NC_000007.14','8':'NC_000008.11','9':'NC_000009.12','10':'NC_000010.11','11':'NC_000011.10','12':'NC_000012.12','13':'NC_000013.11','14':'NC_000014.9','15':'NC_000015.10','16':'NC_000016.10','17':'NC_000017.11','18':'NC_000018.10','19':'NC_000019.10','20':'NC_000020.11','21':'NC_000021.9','22':'NC_000022.11','X':'NC_000023.11','Y':'NC_000024.10','MT':'NC_012920.1'}

RefSeq_accession2Chr = {'NC_000001.11':1,'NC_000002.12':2,'NC_000003.12':3,'NC_000004.12':4,'NC_000005.10':5,'NC_000006.12':6,'NC_000007.14':7,'NC_000008.11':8,'NC_000009.12':9,'NC_000010.11':10,'NC_000011.10':11,'NC_000012.12':12,'NC_000013.11':13,'NC_000014.9':14,'NC_000015.10':15,'NC_000016.10':16,'NC_000017.11':17,'NC_000018.10':18,'NC_000019.10':19,'NC_000020.11':20,'NC_000021.9':21,'NC_000022.11':22,'NC_000023.11':'X','NC_000024.10':'Y','NC_012920.1':'MT'}


def locus2rsID(locusStr,ref=None,alt=None):

    chrN,tok = locusStr.split(':')
    posSta,posEnd = tok.split('-')

    tb_rsID = tabix.open(GCF_path)
    result_records = tb_rsID.querys('%s:%s-%s' % (Chr2RefSeq_accession[chrN],posSta,posEnd))

    resultL = []

    for result in result_records:
        resultL.append((result[:]))
        
    rsIDL = []
    if (ref != None) and (alt != None):
        for i in resultL:
            if (posSta == i[1]) and (ref == i[3]) and (alt in i[4]):
                # rsID = i[2] 
                return i[2]

    else:
        for i in resultL:
            ref = i[3]
            if len(i[4].split(',')) == 1:
                alt = i[4]
                rsID = i[2]
                if len(ref) >= len(alt):
                    locus = ('chr%s:%s-%s' % (chrN,i[1],(int(i[1])+len(ref)-1)))
                elif len(ref) < len(alt):
                    locus = ('chr%s:%s-%s' % (chrN,i[1],(int(i[1])+len(ref)-1)))
                else:
                    raise ValueError
                rsIDL.append((locus, ref, alt, rsID))
            
            else:
                for j in i[4].split(','):
                    alt = j
                    rsID = i[2]
                    if len(ref) >= len(alt):
                        locus = ('chr%s:%s-%s' % (chrN,i[1],(int(i[1])+len(ref)-1)))
                    elif len(ref) < len(alt):
                        locus = ('chr%s:%s-%s' % (chrN,i[1],(int(i[1])+len(ref)-1)))
                    else:
                        raise ValueError
                    rsIDL.append((locus, ref, alt, rsID))
        
        return rsIDL


def rsID2locus(rsID):
    column1 = int(rsID[-1])
    if len(rsID) == 3:
        column2 = 1
        column3 = 1
    else:
        column2 = int(rsID[2:-1])
        column3 = column2
    tb_rsID = tabix.open(GCF_H_path)
    result_records = tb_rsID.querys('%s:%s-%s' % (column1,column2,column3))

    resultL = []

    for result in result_records:
        resultL.append((result[:]))
        
    locusL = []
    for i in resultL:
        ref = i[3]
        alt = i[4]
        locus = i[5]
        locusL.append((locus, ref, alt))
    
    return locusL

