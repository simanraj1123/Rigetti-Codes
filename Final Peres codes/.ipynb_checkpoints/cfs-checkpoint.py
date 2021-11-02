import multiprocessing as mp

def func(x, q):
    print(f'{x**2} from process {x}')
    q.put(x)
    
qout = mp.Queue()

if __name__ == '__main__':
    procs = [mp.Process(target=func, args=(i, qout)) for i in range(10)]
    
    for proc in procs:
        proc.start()
        proc.join()
        
    res = []
    for proc in procs:
        res.append(qout.get())
    print(res)