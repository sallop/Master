#!/usr/bin/python env

def all_permutations(list, length=0):
	def _impl(length, list, ret=[], last=[]):
		if len(last) == length:
			ret.append(last[::])
		else:
			for item in list:
				last.append(item)
				tmpList = list[::]
				tmpList.remove(item)
				_impl(length, tmpList, ret, last)
				last.remove(item)
	ret=[]
	if length is 0 or length > len(list):
		length = len(list)
	_impl(length,list,ret)
	return ret

list = all_permutations([0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.4, 0.1])

for i, term in enumerate(list):
	print "%d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n"%(i, term[0],term[1],term[2],term[3],term[4],term[5],term[6],term[7])

