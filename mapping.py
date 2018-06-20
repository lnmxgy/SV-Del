# _*_ coding:utf-8 _*_
#将第一种情况的read与ref'比对
#命令 ：python mapping.py SAMPLE1.delresult ns150.fq 10 ./result1/ 10（两个10依次为seed大小，候选位置seed支持度）
from collections import defaultdict
from collections import OrderedDict
import re
import os

def make_seed_table(seq, seed_size):
    seed_table = defaultdict(list)#字典键为seed序列，值为该序列出现过的所有位置，为一个列表
    for i in range(0, len(seq) - seed_size):
        seed_table[seq[i:i + seed_size]].append(i)
    return seed_table


def split_read(read, seed_size):
    if len(read) % seed_size:
        msg = "read size:%s can't divided by seed_size:%s." % \
              (len(read), seed_size)
        raise NameError(msg)
    pn = '.{%s}' % seed_size#正则表达式，匹配任意字符，匹配长度为seed大小
    return re.findall(pn, read)#从read中找所有满足正则表达式的序列，返回一个列表，列表中每个元素为一个10长序列（seed）


def candidate(read, seed_table, seed_size):
    seeds = split_read(read, seed_size)
    tmp = {}
    index = 0#表示当前是第几个seed，从0开始
    for seed in seeds:
        if seed_table.has_key(seed):
            for loc in seed_table[seed]:
                loc -= seed_size * index#保证loc为该read的起始位置
                tmp.setdefault(loc, 0) #loc不存在时设置新值
                tmp[loc] += 1  #对一个正确匹配的位置seed支持数进行记录
        index += 1#就算不匹配也要将index加1
    return tmp

def seed_distance(a, b):#计算两个等长序列的距离
    if len(a) != len(b):
        return max(len(a), len(b))
    i = 0
    dif = 0
    while i < len(a):
        if a[i] != b[i]:
            dif += 1
        i += 1
    return dif

def reads_mapping(reads_seq, seed_table, seed_size, h, resultpath, refk, refseq, reads_seq_res, mapdis):#所有read与每一个ref'进行匹配
    readlength = len(reads_seq.values()[0])
    if not os.path.exists(resultpath):
        os.mkdir(resultpath)#创建结果目录
    rst_file = resultpath + refk
    os.mknod(rst_file)#创建结果文件
    f = file(rst_file, 'w')
    reacord = ['read_no', '    |', 'map_information', '    |', 'read_seq', '    |', 'maplocation_ham']
    f.writelines(reacord)
    f.write('\n')
    for k, v in reads_seq.iteritems():
        f.write(k)
        vreverse = ''.join(["ATCG"["TAGC".index(n)] for n in v[::-1]])
        tmp1 = candidate(v, seed_table, seed_size)#用于存储正向比对所有候选位置的支持度
        tmp2 = candidate(vreverse, seed_table, seed_size)#用于存储反向比对所有候选位置的支持度
        candidate_locations1 = [c for c in tmp1 if tmp1[c] > h]
        candidate_locations2 = [c for c in tmp2 if tmp2[c] > h]
        max1 = 0
        max2 = 0
        if len(tmp1) > 0:
            max1 = max(zip(tmp1.values(), tmp1.keys()))[0]#求正向最大支持度
        if len(tmp2) > 0:
            max2 = max(zip(tmp2.values(), tmp2.keys()))[0]#求反向最大支持度
        if max1 > max2:#将方向支持度写入文件
            f.write('    |')
            f.write('z' + str(max1))
        else:
            f.write('    |')
            f.write('r' + str(max2))
        f.write('    |')
        f.write(v)
        f.write('    |')
        if candidate_locations1:
            edisref1loc1list = dict()
            for candidate_location in candidate_locations1:
                if candidate_location > len(refseq) - readlength or candidate_location < 0:#要求read全长序列都在ref'上
                    continue
                edisref1loc1 = seed_distance(v, refseq[candidate_location  : candidate_location + readlength])
                if edisref1loc1 > 60:#计算read与候选位置后150长序列的汉明距离，对候选位置进行验证。
                    continue
                edisref1loc1list[candidate_location] = edisref1loc1
            if len(edisref1loc1list) > 0:
                if mapdis != 0 and edisref1loc1list[edisref1loc1list.keys()[0]] < mapdis:
                    if reads_seq_res.has_key(k[1:]):
                        reads_seq_res[k[1:]] += 1
                    else:
                        reads_seq_res[k[1:]] = 1
                edisref1loc1list = OrderedDict(sorted(edisref1loc1list.items(), key=lambda e:e[1], reverse=False))#字典按value值升序
                f.write(str(edisref1loc1list.keys()[0]))#选汉明距离最小的
                f.write('_' + str(edisref1loc1list[edisref1loc1list.keys()[0]]))
                f.write('_ |')
        if candidate_locations2:
            edisref1loc2list = dict()
            for candidate_location in candidate_locations2:
                if candidate_location > len(refseq) - readlength or candidate_location < 0:
                    continue
                edisref1loc2 = seed_distance(vreverse, refseq[candidate_location : candidate_location + readlength])
                if edisref1loc2 > 60:
                    continue
                edisref1loc2list[candidate_location] = edisref1loc2
            edisref1loc2list = OrderedDict(sorted(edisref1loc2list.items(), key=lambda e:e[1], reverse=False))
            if len(edisref1loc2list) > 0:
                if mapdis != 0 and edisref1loc2list[edisref1loc2list.keys()[0]] < mapdis:
                    if reads_seq_res.has_key(k[1:]):
                        reads_seq_res[k[1:]] += 1
                    else:
                        reads_seq_res[k[1:]] = 1
                f.write(str(edisref1loc2list.keys()[0]))
                f.write('_' + str(edisref1loc2list[edisref1loc2list.keys()[0]]))
                f.write('_*|')
        f.write('\n')
    f.close()





