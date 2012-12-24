#!/usr/bin/env python
import re
from pylab import *

def make_ublas_example(m,n):
    def make_rowstr():
        return "(%s)"%(",".join(['%f'%uniform() for j in range(n)]))

    mat = ",".join([make_rowstr() for i in range(m)])
    return "[%d,%d](%s)"%(m,n,mat)

def read_ublas_mat(line):
    """
    [,]((,,),(,,),(,,))
    [rowsz,colsz](rowvecs)
    """
    fmt_mat = re.compile("\[(\d+),(\d+)\]\((.+)\)")
    (rowsz,colsz,rowvecs) = re.compile(fmt_mat).search(line).groups()
    rowsz, colsz = int(rowsz), int(colsz)
    fmt_rows = ",".join(["\((.+)\)"]*rowsz)
    peeled_rows = re.compile(fmt_rows).search(rowvecs).groups()
    rows = map(lambda s: map(lambda x:float(x), s.split(",")), peeled_rows)
    return rows

def read_ublas_vec(line):
    """
    [](,,)
    [size](elems)
    """
    fmt_vec = re.compile("\[(\d+)\]\((.+)\)")
    (vecsz,elems) = re.compile(fmt_vec).search(line).groups()
    vecsz = int(vecsz)
    vec = map(lambda s:float(s),elems.split(","))
    return vec

def read_ublas_matstr(line):
    """
    [,]((,,),(,,),(,,))
    [rowsz,colsz](rowvecs)
    """
    fmt_mat = re.compile("\[(\d+),(\d+)\]\((.+)\)")
    (rowsz,colsz,rowvecs) = re.compile(fmt_mat).search(line).groups()
    rowsz, colsz = int(rowsz), int(colsz)
    fmt_rows = ",".join(["\((.+)\)"]*rowsz)
    peeled_rows = re.compile(fmt_rows).search(rowvecs).groups()
    rows = map(lambda s:s.split(","), peeled_rows)
    return rows

def read_ublas_vecstr(line):
    """
    [](,,)
    [size](elems)
    """
    fmt_vec = re.compile("\[(\d+)\]\((.+)\)")
    (vecsz,elems) = re.compile(fmt_vec).search(line).groups()
    vecsz = int(vecsz)
    vec = elems.split(",")
    return vec


def read_t_mat(line):
    """
    -mat.dat
    t    [,]((,,),(,,),(,,))
    t    [rowsz,colsz](rowvecs)
    return (t,rows)
    """
    re_t = re.compile("^(\d+)").search(line)
    time = int(re_t.group(1))
    (s,e) = re_t.span()
    rows = read_ublas_mat(line[e:])
    return (time,rows)

def read_t_gibbs(line):
    """
    -t_gibbs.dat
    t    [](,,)
    t    [size](elems)
    return (t,vec)
    """
    re_t = re.compile("^(\d+)").search(line)
    time = int(re_t.group(1))
    (s,e) = re_t.span()
    vec = read_t_gibbs(line[e:])
    return (time,vec)

def read_t_matstr(line):
    """
    -mat.dat
    t    [,]((,,),(,,),(,,))
    t    [rowsz,colsz](rowvecs)
    return (t,rows)
    """
    re_t = re.compile("^(\d+)").search(line)
    time = re_t.group(1)
    (s,e) = re_t.span()
    rows = read_ublas_mat(line[e:])
    return (time,rows)

def read_t_gibbsstr(line):
    """
    -t_gibbs.dat
    t    [](,,)
    t    [size](elems)
    return (t,vec)
    """
    re_t = re.compile("^(\d+)").search(line)
    time = re_t.group(1)
    (s,e) = re_t.span()
    vec = read_t_gibbs(line[e:])
    return (time,vec)


if __name__ == '__main__':
    fr = open(sys.argv[1],'r')
    for line in fr:
        #read_ublas(line)
        (t,mat) = read_t_mat(line)

        #print "%d\t"%(t),
        print "%s\t"%(t),
        for row in mat:
            for elem in row:
                print "%s\t"%(elem),

        print
    fr.close()
