import sys
from multiprocessing import Process, Queue
from cutadapt.scripts.cutadapt import AdapterCutter
from cutadapt.modifiers import QualityTrimmer, UnconditionalCutter
from cutadapt.seqio import FastqReader
import cutadapt
from distutils.version import StrictVersion

cutadapt_version = StrictVersion(cutadapt.__version__)
if cutadapt_version < StrictVersion('1.8.1'):
	raise ImportError('miRge requires cutadapt >= version 1.8.1 to operate.')
elif cutadapt_version < StrictVersion('1.10.0'):
	from cutadapt.adapters import Adapter, gather_adapters
	def parse_adapters(adapter, error_rate=None):
		adapters = []
		for name,seq,where in gather_adapters(adapter.split(','), [], []):
			adapters.append(Adapter(seq, where, error_rate, name=name))
		return adapters
else:
	from cutadapt.adapters import AdapterParser
	def parse_adapters(adapter, error_rate=None):
		return AdapterParser(max_error_rate=error_rate).parse_multi(adapter.split(','), [], [])

class Worker(Process):
	def __init__(self, queue=None, results=None, adapter=None, phred64=False):
		super(Worker, self).__init__()
		self.queue=queue
		self.results = results
		self.phred = 64 if phred64 else 33
		self.modifiers = [QualityTrimmer(0, 10, self.phred)]
		self.adapters = []
		self.error_rate = 0.12
		self.min_length = 16
		if adapter.startswith('+'):
			self.modifiers.append(UnconditionalCutter(int(adapter)))
		elif adapter == 'none':
			self.adapter = None
		else:
			self.adapters = parse_adapters(adapter, error_rate=self.error_rate)
			adapter_cutter = AdapterCutter(self.adapters)
			self.modifiers.append(adapter_cutter)

	def run(self):
		# we can't use the sentinel iter(self.queue.get, None) because of some issue with Sequence classes
		results = self.results
		get_func = self.queue.get
		reads = get_func()
		modifiers = self.modifiers
		min_length = self.min_length
		while reads is not None:
			result_batch = []
			for read in reads:
				for modifier in modifiers:
					read = modifier(read)
				if len(read.sequence) >= min_length:
					try:
						second_header = read.name2
					except AttributeError:
						second_header = read.name if getattr(read, 'twoheaders', False) else ''
					if second_header:
						read_out = '@' + read.name + '\n' + read.sequence + '\n+' + second_header + '\n' +  read.qualities + '\n'
					else:
						read_out = '@' + read.name + '\n' + read.sequence + '\n+\n' + read.qualities + '\n'
					result_batch.append(read_out)
			results.put(result_batch)
			reads = get_func()

class Writer(Process):
	def __init__(self, queue=None, trimmed=None, outfile=None):
		super(Writer, self).__init__()
		self.queue = queue
		self.trimmed = trimmed
		self.outfile = outfile

	def run(self):
		get_func = self.queue.get
		reads = get_func()
		with open(self.outfile, 'w') as outfile:
		#outfile = self.outfile
			kept = 0
			while reads is not None:
				for read in reads:
					outfile.write(read)
					kept += 1
				reads = get_func()
		self.trimmed.put(kept)


def trim_file(infile, adapter, outfile, threads=1, phred=33):
	read_queue = Queue()
	result_queue = Queue()
	trimmed_queue = Queue()
	workers = []
	def start_workers():
		for i in xrange(threads):
			worker = Worker(queue=read_queue, results=result_queue, phred64=phred==64, adapter=adapter)
			workers.append(worker)
			worker.start()
	writer = Writer(queue=result_queue, trimmed=trimmed_queue, outfile=outfile)
	writer.start()
	batch = []
	for index, read in enumerate(FastqReader(infile)):
		batch.append(read)
		if index < 1000 and phred == 33:
			if any([i for i in read.qualities if ord(i) > 74]):
				phred = 64
		if index % 10000 == 0:
			if not workers:
				start_workers()
			read_queue.put(batch)
			batch = []
	if not workers:
		start_workers()
	read_queue.put(batch)
	processed = index+1
	# poison pill to stop workers
	for i in xrange(threads):
		read_queue.put(None)

	for i in workers:
		i.join()

	# poison pill for writers
	result_queue.put(None)

	# wait for writing to finish
	writer.join()
	#print "Output done"

	trimmed_queue.put(None)

	kept_reads = sum([i for i in iter(trimmed_queue.get, None)])
	
	return (phred, processed, kept_reads)
	#with logfile as o:
	#	o.write('Starting reads: {0}\n'.format(processed))
	#	o.write('Processed reads: {0}\n'.format(kept_reads))
	#print ('{0}\n'.format(phred))