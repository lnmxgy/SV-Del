# _*_ coding:utf-8 _*_
#识别deletion

from mapping import *

def readresult(simvaluefile, listBoth):#读取比对后的信息进行处理
    with open(simvaluefile, 'r') as fp:
        fp.readline()
        while True:
            line1 = fp.readline().strip('\n')
            line2 = fp.readline().strip('\n')
            if not line1:
                break
            if not line2:
                break
            tmp1 = line1.split("|")
            tmp2 = line2.split("|")
            temp1 = tmp1[3].split("_")
            temp2 = tmp2[3].split("_")
            if len(tmp1) > 4 and len(tmp2) > 4:
                if int(temp1[1]) >2 or int(temp2[1]) > 2:
                    if int(temp1[1]) < int(temp2[1]):
                        listi = [tmp2[0], tmp2[2], tmp2[3]]
                        listBoth.append(listi)
                    else:
                        listi = [tmp1[0], tmp1[2], tmp1[3]]
                        listBoth.append(listi)


def buildbreakup(start, end, breakp, ref, readname):#根据断点长度+位置创建哈希表，长度不为0
    left = start
    right = start
    if ref[start] == ref[end]:
        i1 = 0
        while(end + i1 < len(ref) and ref[start + i1] == ref[end + i1]):
            i1 += 1
        right = start + i1
    else:
        i2 = 1
        while(start - i2 >= 0 and ref[start - i2] == ref[end - i2]):
            i2 += 1
        left = start - i2 + 1
    if left == right:
        key = str(start) + " " + str(end - start)
    else:
        key = str(left) + " " + str(right) + " " + str(end - start)
    breakp[key].append(readname)



def findLongestCommonSubtringWithTolerance(s1, s2, tolerance):
    dp = [[[0 for j in range(len(s2)+1)] for i in range(len(s1)+1)]  for k in range(tolerance + 1)]

    fmax = 0

    for k in range(tolerance + 1):
        for i in range(1, len(s1) + 1):
            for j in range(1, len(s2) + 1):
                if s1[i - 1] == s2[j - 1]:
                    dp[k][i][j] = dp[k][i - 1][j - 1] + 1
                else:
                    if k > 0:
                        dp[k][i][j] = dp[k-1][i - 1][j - 1] + 1
                    else:
                        dp[k][i][j] = 0
                if dp[k][i][j] > fmax:
                    fmax = dp[k][i][j]

    candidate = []
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            if dp[tolerance][i][j] == fmax:
                candidate.append([fmax, i - fmax, j - fmax])
    return candidate


def findBp(read, ref, readname, breakp, readlength, tolerance, mapl):
    reff = ref[mapl : mapl + readlength]
    i = 0
    lmmax = 0
    toleranceused = 0
    canl = []
    while(i <= tolerance):
        candidatel = findLongestCommonSubtringWithTolerance(read, reff, i)
        mmax = candidatel[0][0]
        if mmax - lmmax > 5:
            if lmmax != 0:
                toleranceused += 1
            lmmax = mmax
            canl = candidatel
        else:
            break
        i += 1
    if len(canl) == 1:
        lmmax, lreadstart, lrefstart = canl[0]
        if lmmax == len(read):
            return
        lrefstart += mapl
        if lreadstart == 0:
            refleft = ref[lrefstart + lmmax: ]
            readleft = read[lmmax : ]
            if len(readleft) < 10:
                return
            i = 0
            canr = []
            rmmax = 0
            while(i <= tolerance - toleranceused):
                candidater = findLongestCommonSubtringWithTolerance(readleft, refleft, i)
                mmax = candidater[0][0]
                if mmax - rmmax > 5:
                    canr = candidater
                    rmmax = mmax
                else:
                    break
                i += 1
            if rmmax != len(readleft):
                return
            for canri in canr:
                rmmax, rreadstart, rrefstart = canri
                rrefstart = rrefstart + lrefstart + lmmax
                start = lrefstart + lmmax
                end = rrefstart
                buildbreakup(start, end, breakp, ref, readname)
        elif lreadstart + lmmax == len(read):
            refleft = ref[0 : lrefstart]
            readleft = read[0 : readlength - lmmax]
            if len(readleft) < 10:
                return
            i = 0
            canr = []
            rmmax = 0
            while(i <= tolerance - toleranceused):
                candidater = findLongestCommonSubtringWithTolerance(readleft, refleft, i)
                mmax = candidater[0][0]
                if mmax - rmmax > 5:
                    canr = candidater
                    rmmax = mmax
                else:
                    break
                i += 1
            if rmmax != len(readleft):
                return
            for canri in canr:
                rmmax, rreadstart, rrefstart = canri
                start = rrefstart + rmmax
                end = lrefstart
                buildbreakup(start, end, breakp, ref, readname)
        else:
            print "in"


def findBreak(refSimvalueFilePath, k, v, u, tolerance):
    print k
    breakp = defaultdict(list)#断点信息
    listBoth = []#存两端都有比对位置，汉明距离大的一端的信息
    readresult(refSimvalueFilePath + k, listBoth)
    if(len(listBoth) == 0):
        return breakp, 0
    ref = v
    halfreflen = (len(ref) - (int(k.split(".")[0].split("_")[1]))) / 2
    readlength = len(listBoth[0][1].strip())
    print len(listBoth)
    for i in range(0, len(listBoth)):#根据listBoth信息预测断点
        mapl = int(listBoth[i][2].split("_")[0])
        if '*' in listBoth[i][2]:#read反向匹配到ref'上
            read = ''.join(["ATCG"["TAGC".index(n)] for n in listBoth[i][1].strip()[::-1]])
        else:
            read = listBoth[i][1].strip()
        findBp(read, ref, listBoth[i][0].strip(), breakp, readlength, tolerance, mapl)


    bpk = list(breakp.keys())
    for j in range(0, len(breakp)):
        for l in range(j + 1, len(breakp)):
            breakp[bpk[l]] = [m for m in breakp[bpk[l]] if m not in breakp[bpk[j]]]#使一个read只支持一个断点

    breakpB = OrderedDict(sorted(breakp.items(), key=lambda e:len(e[1]), reverse=True))#对一个断点得到的deletion进行排序
    beakpdel = []
    for key in bpk:
        if len(breakpB[key]) < 2:#删除支持度过小的断点
            beakpdel.append(key)
    for key in beakpdel:
        del breakpB[key]

    for ii, jj in breakpB.iteritems():
        print ii, len(jj)

    return breakpB, 1


def selectBreak(nonemapref, dictref, breakpfe, num):
    for nonkey in nonemapref:
        del dictref[nonkey]
    breakpf = OrderedDict()
    for i in range(0, len(breakpfe)):
        k = dictref.keys()[i]
        v = dictref[k]
        counttop4 = 0
        halfreflen = (len(v) - (int(k.split(".")[0].split("_")[1]))) / 2
        for key, value in breakpfe[i].iteritems():#构建deletion真实位置的字典，最多取k个
            if counttop4 == num:
                break
            tmp = key.split()
            if len(tmp) > 2:
                key1 = int(tmp[0]) + int(k.split("_")[0]) - halfreflen
                key2 = int(tmp[1]) + int(k.split("_")[0]) - halfreflen
                keyf = str(key1) + "_" + str(key2) + " " + tmp[2]
                if breakpf.has_key(keyf):
                    print "close deletion exist"
                    print keyf
                else:
                    breakpf[keyf] = value
            else:
                keyf = str(int(tmp[0]) + int(k.split("_")[0]) - halfreflen) + " " + tmp[1]
                if breakpf.has_key(keyf):
                    print "close deletion exist"
                else:
                    breakpf[keyf] = value
            counttop4 += 1
    writecase1result(breakpf)


def writecase1result(breakpf):
    print "case1write"
    fp = open("BreakPoint.txt", 'w')
    fp.write('startPos(or start range)    delent    supports\n')
    for k, v in breakpf.iteritems():
        fp.write(k)
        fp.write(" " + str(len(v)))
        fp.write('\n')
    fp.close()

def writeCase2ReadFile(dictref, nonemapref, refSimvalueFilePath, ns150):#生成第二种read
    readnum = []
    for k in dictref.keys():#得到与第一种情况比对的read号
        if k not in nonemapref:
            with open(refSimvalueFilePath + k, 'r') as fp:
                fp.readline()
                while True:
                    line1 = fp.readline().strip('\n')
                    line2 = fp.readline().strip('\n')
                    if not line1:
                        break
                    if not line2:
                        break
                    tmp1 = line1.split("|")
                    tmp2 = line2.split("|")
                    if len(tmp1) > 4 or len(tmp2) > 4:
                        readnum.append(tmp1[0].strip())
                        readnum.append(tmp2[0].strip())
            fp.close()
    for i in range(0, len(readnum)):
        if ns150.has_key(readnum[i]):#ns150中删除第一种情况的read
            del ns150[readnum[i]]
    return ns150

def gecase2ref(oref, case2bpInfo, caseIIbpInfo, insertsize, readlength, lenl):#根据参考基因生成第二种类型的ref'，
    case2ref = OrderedDict()
    for i in range(0, len(case2bpInfo)):
        tmp = case2bpInfo[i].split(".")[0].split("_")
        refk = str(tmp[0]) + '_' + str(tmp[1]) + ".smvalue"
        if caseIIbpInfo[case2bpInfo[i]] == "none":
            startp = int(tmp[0])
            delen = int(tmp[1])
        else:
            tp = caseIIbpInfo[case2bpInfo[i]].split()
            startp = int(tp[0])
            delen = int(tp[1])
        nref = oref[startp - lenl : startp] + oref[startp + delen : startp + delen + lenl]
        case2ref[refk] = nref
    return case2ref




