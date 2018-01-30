import os
from reportlab.graphics.charts.textlabels import Label
from reportlab.graphics.shapes import Drawing, _DrawingEditorMixin 
from reportlab.graphics.charts.barcharts import VerticalBarChart, HorizontalBarChart
from reportlab.lib.units import inch 
from reportlab.platypus.flowables import Image, Spacer, PageBreak
from reportlab.platypus.paragraph import Paragraph 
from reportlab.platypus.doctemplate import SimpleDocTemplate 
from reportlab.lib.styles import getSampleStyleSheet 
from reportlab.platypus.tables import Table, TableStyle, GRID_STYLE, BOX_STYLE, LABELED_GRID_STYLE, COLORED_GRID_STYLE, LIST_STYLE, LongTable 
from reportlab.lib import colors
from reportlab.lib.formatters import DecimalFormatter

class TableBarChart1(_DrawingEditorMixin,Drawing): 
	def __init__(self,xList,yList,width=400,height=200, *args,**kw): 
		Drawing.__init__(self,width,height, *args,**kw) 
		self.width = 80*1.2*2.5
		self.height = 60*0.8*2.5
		self._add(self,HorizontalBarChart(),name='chart',validate=None,desc=None) 
		self.chart.y = 11.2*2
		self.chart.x = 20*2
		#self.chart.data = [[0.05, 0.05, 0.1, 0.1, 0.05, 0.05, 0.955555]]
		self.chart.data = yList
		self.chart.width = self.width*0.75
		self.chart.height = self.height*0.75
		self.chart.bars.strokeColor = None
		self.chart.bars[0].fillColor = colors.blue
		self.chart.valueAxis.valueMin = 0
		self.chart.valueAxis.valueMax = 1.0
		self.chart.valueAxis.valueStep = 0.2
		self.chart.categoryAxis.visibleTicks = False
		self.chart.barLabelFormat = DecimalFormatter(3)
		self.chart.barLabels.boxAnchor = 'w'
		self.chart.barLabels.fontSize = 2.8*2.5
		self.chart.barLabels.fontName = 'Helvetica-Bold'
		self.chart.valueAxis.labels.fontSize = 2.8*2.5
		self.chart.valueAxis.labels.fontName = 'Helvetica-Bold'
		#self.chart.categoryAxis.categoryNames = ['unaligned', 'mRNA', 'rRNA', 'snoRNA', 'tRNA', 'hairpin miRNA', 'miRNA']
		self.chart.categoryAxis.categoryNames = xList
		self.chart.categoryAxis.labels.fontSize = 2.8*2.5
		self.chart.categoryAxis.labels.fontName = 'Helvetica-Bold'
		self._add(self,Label(),name='XLabel',validate=None,desc='Percentage')
		self.XLabel.fontSize = 4.4*2
		self.XLabel.fontName = 'Helvetica-Bold'
		self.XLabel._text = 'Percentage'
		self.XLabel.dx = 48*2.5
		self.XLabel.dy = 2

class TableBarChart2(_DrawingEditorMixin,Drawing):
	def __init__(self,xList,yList,width=400,height=200, *args,**kw):
		Drawing.__init__(self,width,height,*args,**kw) 
		self.width = 80*1.2*2.5
		self.height = 60*0.8*2.5
		self._add(self,VerticalBarChart(),name='chart',validate=None,desc=None) 
		self.chart.y = 11.2*2
		self.chart.x = 20*2
		#self.chart.data = [[0.05, 0.05, 0.1, 0.1, 0.05, 0.05, 0.955555]]
		self.chart.data = yList
		self.chart.width = self.width*0.75
		self.chart.height = self.height*0.75
		self.chart.bars.strokeColor = None
		self.chart.bars[0].fillColor = colors.blue
		maxTmp = (max(yList[0])/1000)*1000
		valueStepsTmp = []
		for i in [0.2, 0.4, 0.6, 0.8, 1.0]:
			valueStepsTmp.append(i*maxTmp)
		self.chart.valueAxis.valueSteps = valueStepsTmp
		self.chart.categoryAxis.visibleTicks = False
		self.chart.barLabels.fontSize = 2.8*2.5
		self.chart.barLabels.fontName = 'Helvetica-Bold'
		self.chart.valueAxis.labels.fontSize = 2.8*2.5
		self.chart.valueAxis.labels.fontName = 'Helvetica-Bold'
		xListNew = []
		interval = (max([int(item) for item in xList])-min([int(item) for item in xList]))/10
		if interval == 0:
			interval = (max([int(item) for item in xList])-min([int(item) for item in xList]))/5
			if interval == 0:
				interval = max([int(item) for item in xList])-min([int(item) for item in xList])
		for i in range(len(xList)):
			if i%interval == 0:
				xListNew.append(xList[i])
			else:
				xListNew.append('')
		self.chart.categoryAxis.categoryNames = xListNew
		self.chart.categoryAxis.labels.fontSize = 2.8*2.5
		self.chart.categoryAxis.labels.fontName = 'Helvetica-Bold'
		self._add(self,Label(),name='XLabel',validate=None,desc='Percentage')
		self.XLabel.fontSize = 4.4*2
		self.XLabel.fontName = 'Helvetica-Bold'
		self.XLabel._text = 'Read length'
		self.XLabel.dx = 48*2.5
		self.XLabel.dy = 2
		self._add(self,Label(),name='YLabel',validate=None,desc='Percentage')
		self.YLabel.fontSize = 4.4*2
		self.YLabel.fontName = 'Helvetica-Bold'
		self.YLabel._text = 'Counts'
		self.YLabel.dx = 0.8
		self.YLabel.dy = 32*2.5
		self.YLabel.angle = 90

def generateReport(outputdir, sampleList, readLengthDic, logDic, annotNameList, seqDic, spikeIn):
	filename = os.path.join(outputdir, 'annotationReport.pdf')
	annotationFile = os.path.join(outputdir, 'annotation.report.csv')
	outf = open(annotationFile, 'w')
	outf.write('File name(s),Total Input Reads,Trimmed Reads(all),Trimmed Reads(unique),All miRNA Reads,Filtered miRNA Reads,Unique miRNAs,Hairpin miRNAs,tRNA Reads,snoRNA Reads,rRNA Reads,ncRNA others,mRNA Reads,Remaining Reads\n')
	styleSheet = getSampleStyleSheet()
	story = []
	if spikeIn:
		data1_1= [[Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>File name(s)</font></b></para>",styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Total Input Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Trimmed Reads(all / unique)</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>All miRNA Reads / Filtered miRNA Reads</font></b></para>",styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Unique miRNAs</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Hairpin miRNAs</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>tRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>snoRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>rRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>ncRNA others Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>mRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>spike-in Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Remaining Reads</font></b></para>", styleSheet["BodyText"])]]
	else:
		data1_1= [[Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>File name(s)</font></b></para>",styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Total Input Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Trimmed Reads(all / unique)</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>All miRNA Reads / Filtered miRNA Reads</font></b></para>",styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Unique miRNAs</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Hairpin miRNAs</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>tRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>snoRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>rRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>ncRNA others Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>mRNA Reads</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Remaining Reads</font></b></para>", styleSheet["BodyText"])]]
	data1_2 = [[Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>File name(s)</font></b></para>",styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Reads composition distrubution</font></b></para>", styleSheet["BodyText"]),
				Paragraph("<para align=center><b><font color=black size=12.0 face='Helvetica-Bold'>Reads length distrubution</font></b></para>", styleSheet["BodyText"])]]
	for i in range(len(sampleList)):
		if spikeIn:
			xList1 = ['unaligned', 'mRNA', 'ncRNA others', 'rRNA', 'snoRNA', 'tRNA', 'hairpin miRNA', 'miRNA', 'spike-in']
		else:
			xList1 = ['unaligned', 'mRNA', 'ncRNA others', 'rRNA', 'snoRNA', 'tRNA', 'hairpin miRNA', 'miRNA']
		yList1 = [[]]
		yList1[0].append(float(logDic['quantStats'][i]['remReads'])/logDic['quantStats'][i]['trimmedReads'])
		yList1[0].append(float(logDic['quantStats'][i]['mrnaReads'])/logDic['quantStats'][i]['trimmedReads'])
		yList1[0].append(float(logDic['quantStats'][i]['ncrnaOthersReads'])/logDic['quantStats'][i]['trimmedReads'])
		yList1[0].append(float(logDic['quantStats'][i]['rrnaReads'])/logDic['quantStats'][i]['trimmedReads'])
		yList1[0].append(float(logDic['quantStats'][i]['snornaReads'])/logDic['quantStats'][i]['trimmedReads'])
		yList1[0].append(float(logDic['quantStats'][i]['trnaReads'])/logDic['quantStats'][i]['trimmedReads'])
		yList1[0].append(float(logDic['quantStats'][i]['hairpinReads'])/logDic['quantStats'][i]['trimmedReads'])
		yList1[0].append(float(logDic['quantStats'][i]['mirnaReads'])/logDic['quantStats'][i]['trimmedReads'])
		if spikeIn:
			yList1[0].append(float(logDic['quantStats'][i]['spikeInReads'])/logDic['quantStats'][i]['trimmedReads'])
		B1 = TableBarChart1(xList1,yList1) 
		xList2 = []
		yList2 = [[]]
		tmpList = readLengthDic.keys()
		tmpList.sort()
		for readLen in tmpList:
			xList2.append(str(readLen))
			yList2[0].append(readLengthDic[readLen][i])
		B2 = TableBarChart2(xList2,yList2)
		dataTmp1_1 = []
		dataTmp1_1.append(sampleList[i])
		dataTmp1_1.append(str(logDic['quantStats'][i]['totalReads']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['trimmedReads'])+' / '+str(logDic['quantStats'][i]['trimmedUniq']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['mirnaReads'])+' / '+str(logDic['quantStats'][i]['mirnaReadsFiltered']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['mirnaUniqFiltered']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['hairpinReads']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['trnaReads']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['snornaReads']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['rrnaReads']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['ncrnaOthersReads']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['mrnaReads']))
		if spikeIn:
			dataTmp1_1.append(str(logDic['quantStats'][i]['spikeInReads']))
		dataTmp1_1.append(str(logDic['quantStats'][i]['remReads']))
		dataTmp1_2 = []
		dataTmp1_2.append(sampleList[i])
		dataTmp1_2.append(B1)
		dataTmp1_2.append(B2)
		data1_1.append(dataTmp1_1)
		data1_2.append(dataTmp1_2)
		outf.write(sampleList[i]+',')
		if spikeIn:
			outf.write(','.join([str(logDic['quantStats'][i]['totalReads']), str(logDic['quantStats'][i]['trimmedReads']), str(logDic['quantStats'][i]['trimmedUniq']), str(logDic['quantStats'][i]['mirnaReads']),
					str(logDic['quantStats'][i]['mirnaReadsFiltered']), str(logDic['quantStats'][i]['mirnaUniqFiltered']), str(logDic['quantStats'][i]['hairpinReads']), str(logDic['quantStats'][i]['trnaReads']),
					str(logDic['quantStats'][i]['snornaReads']), str(logDic['quantStats'][i]['rrnaReads']), str(logDic['quantStats'][i]['ncrnaOthersReads']), str(logDic['quantStats'][i]['mrnaReads']),  str(logDic['quantStats'][i]['spikeInReads']), str(logDic['quantStats'][i]['remReads'])]))
		else:
			outf.write(','.join([str(logDic['quantStats'][i]['totalReads']), str(logDic['quantStats'][i]['trimmedReads']), str(logDic['quantStats'][i]['trimmedUniq']), str(logDic['quantStats'][i]['mirnaReads']),
					str(logDic['quantStats'][i]['mirnaReadsFiltered']), str(logDic['quantStats'][i]['mirnaUniqFiltered']), str(logDic['quantStats'][i]['hairpinReads']), str(logDic['quantStats'][i]['trnaReads']),
					str(logDic['quantStats'][i]['snornaReads']), str(logDic['quantStats'][i]['rrnaReads']), str(logDic['quantStats'][i]['ncrnaOthersReads']), str(logDic['quantStats'][i]['mrnaReads']), str(logDic['quantStats'][i]['remReads'])]))
		outf.write('\n')
	outf.close()
	if spikeIn:
		colCount1_1 = 14-2
		colWidthsArgu = (1.367*inch, None, None, None, None, None, None, None, None, None, None, None, None)
	else:
		colCount1_1 = 13-2
		colWidthsArgu = (1.367*inch, None, None, None, None, None, None, None, None, None, None, None)
	t1=Table(data1_1, colWidths=colWidthsArgu,
		style=[
		('FONTSIZE', (0, 0), (-1, -1), 9),
		('TEXTFONT', (0, 0), (-1, -1), 'Helvetica-Bold'),
		('BACKGROUND', (0, 0), (colCount1_1, 0), colors.lightblue),
		('BOX',(0,0),(-1,-1),1,colors.black),
		('GRID',(0,0),(-1,-1),0.5,colors.black),
		('VALIGN',(0,0),(-1,-1),'MIDDLE'),
		('ALIGN',(0,0),(-1,-1),'CENTER'),
		('LEFTPADDING', (0,0),(-1,-1),4),
		('RIGHTPADDING', (0,0),(-1,-1),4),
		('BOTTOMPADDING', (0,0),(-1,-1),6),
		('TOPPADDING', (0,0),(-1,-1),6),
		])
	colCount1_2 = 2
	t2=Table(data1_2, colWidths=(1.367*inch, None, None),
		style=[
		('FONTSIZE', (0, 0), (-1, -1), 9),
		('TEXTFONT', (0, 0), (-1, -1), 'Helvetica-Bold'),
		('BACKGROUND', (0, 0), (colCount1_2, 0), colors.lightblue),
		('BOX',(0,0),(-1,-1),1,colors.black),
		('GRID',(0,0),(-1,-1),0.5,colors.black),
		('VALIGN',(0,0),(-1,-1),'MIDDLE'),
		('ALIGN',(0,0),(-1,-1),'CENTER'),
		('LEFTPADDING', (0,0),(-1,-1),4),
		('RIGHTPADDING', (0,0),(-1,-1),4),
		('BOTTOMPADDING', (0,0),(-1,-1),6),
		('TOPPADDING', (0,0),(-1,-1),6),
		])
	story.append(Paragraph("<para align=left><b><font color=darkgray size=20 face='Helvetica-Bold'>Sample Result(s)</font></b></para>", styleSheet['Title']))
	story.append(t1)
	story.append(Paragraph("<para align=left><b><font color=darkgray size=20 face='Helvetica-Bold'><br /><br /><br />Read Length and Composition Figures</font></b></para>", styleSheet['Title']))
	story.append(t2)
	story.append(PageBreak())
	story.append(Paragraph("<para align=left><b><font color=darkgray size=20 face='Helvetica-Bold'>Annotation summary of unique sequences from sample set and processing time</font></b></para>", styleSheet['Title']))
	data2 = [[Paragraph("<para align=center spaceb=3><b><font color=black size=12.0 face='Helvetica-Bold'>Annotation-Round</font></b></para>",styleSheet["BodyText"]),
			Paragraph("<para align=center spaceb=3><b><font color=black size=12.0 face='Helvetica-Bold'>Unique Seqs</font></b></para>", styleSheet["BodyText"]),
			Paragraph("<para align=center spaceb=3><b><font color=black size=12.0 face='Helvetica-Bold'>cpuTime(sec)</font></b></para>", styleSheet["BodyText"])]]
	data2.append(['all sequences', str(len(seqDic.keys())), ''])
	
	dataTmp = []
	dataTmp.append(annotNameList[0])
	dataTmp.append(str(logDic['annotStats'][0]['readsAligned']))
	dataTmp.append('%.2f'%(logDic['annotStats'][0]['cpuTime']))
	data2.append(dataTmp)
	dataTmp = []
	dataTmp.append(annotNameList[len(annotNameList)-1])
	dataTmp.append(str(logDic['annotStats'][len(annotNameList)-1]['readsAligned']))
	dataTmp.append('%.2f'%(logDic['annotStats'][len(annotNameList)-1]['cpuTime']))
	data2.append(dataTmp)
	for i in range(1, len(annotNameList)-1):
		dataTmp = []
		dataTmp.append(annotNameList[i])
		dataTmp.append(str(logDic['annotStats'][i]['readsAligned']))
		dataTmp.append('%.2f'%(logDic['annotStats'][i]['cpuTime']))
		data2.append(dataTmp)
	t2=Table(data2,style=[
		('FONTSIZE', (0, 0), (-1, -1), 9),
		('TEXTFONT', (0, 0), (-1, -1), 'Helvetica-Bold'),
		('BACKGROUND', (0, 0), (2, 0), colors.lightblue),
		('BOX',(0,0),(-1,-1),1,colors.black),
		('GRID',(0,0),(-1,-1),0.5,colors.black),
		('VALIGN',(0,0),(-1,-1),'MIDDLE'),
		('ALIGN',(0,0),(-1,-1),'CENTER'),
		('LEFTPADDING', (0,0),(-1,-1),4),
		('RIGHTPADDING', (0,0),(-1,-1),4),
		('BOTTOMPADDING', (0,0),(-1,-1),6),
		('TOPPADDING', (0,0),(-1,-1),6),
		])
	story.append(t2)
	SimpleDocTemplate(filename, pagesize=(984,695.8)).build(story) 

if __name__=='__main__': 
	main() 