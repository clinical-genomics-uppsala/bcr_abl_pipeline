from pysam import VariantFile
import xlsxwriter
from datetime import date
import gzip
import subprocess

vcf = VariantFile(snakemake.input.vcf)

sample = str(vcf.header.samples[0])
today = date.today()
emptyList = ['', '', '', '', '', '']


depth_dict={}
with gzip.open(snakemake.input.mosdepth_regions, 'rt') as mosdepthfile:
    for lline in mosdepthfile:
        line = lline.strip().split("\t")
        depth_dict[line[3]] = int(float(line[4]))

# All positions in Branford list
# Add vaiant change in bedfile chr, pos, mutation name, , ref, alt, c-name
branfordVariants = []
with open(snakemake.input.branford, 'r') as branfordfile:
    for line in branfordfile:
        branfordVariants.append(line.strip().split("\t"))

branfordSNV_dup = []
branfordHeader = ["Mutation", "Transcript", "Chr", "Pos", "Ref", "Alt", "Allelic Freq", "Allelic Depth", "Depth"]

allSNVs = []
allHeader = ["Sample", "Chr", "Pos", "Ref", "Alts", "Allelic Freq", "Depth", "Allelic Depth", "Filter flag"]

for record in vcf.fetch():
    if record.filter.keys() == ["PASS"]:
        filter = ''
    else:
        filter = ", ".join(list(record.filter))

    outline = [sample, record.contig, str(record.pos), record.ref, ";".join(record.alts), ';'.join(map(str, record.info["AF"])), str(record.info["DP"]), ";".join(map(str, record.samples[sample].get("AD"))), filter]
    allSNVs.append(outline)

    for b_line in branfordVariants: #chr, pos, ref, alt, mutation, name
        if record.contig == b_line[0] and record.pos == int(b_line[2]) and  b_line[5] in record.alts:
            dp = str(record.info["DP"])
            index = list(record.alts).index(b_line[5])
            af = str(round(record.info["AF"][index], 3))
            a_dp = str(record.samples[sample].get("AD")[index+1])
            outline = [b_line[6], b_line[3], b_line[0], str(b_line[2]), b_line[4], b_line[5], af, a_dp, dp]
            branfordSNV_dup.append(outline)
            break
        else:
            outline = [b_line[6], b_line[3], b_line[0], str(b_line[2]), b_line[4], b_line[5], "", "", ""]
            branfordSNV_dup.append(outline)

# Remove duplicated lines and add dp from mosdepth on non-call regions
branfordSNV = []
for line in branfordSNV_dup:
    if line not in branfordSNV:
        if line[8] == "":
            line[8] = str(depth_dict[line[1]])
        branfordSNV.append(line)


fusion_lines = []
fusionHeading = ["Gene1", "Gene2", "Gene1_breakpoint", "Gene2_breakpoint", "BCR-ABL1 Type", "Fusion Type", "Split reads1", "Split reads2", "Discordant Mates", "Coverage gene1", "Coverage gene2", "Confidence"]
with open(snakemake.input.arriba_tsv, 'r') as tsv_arriba:
    headerIndex=tsv_arriba.readline().strip().split("\t")
    for lline in tsv_arriba:
        line = lline.strip().split("\t")
        if line[headerIndex.index("confidence")] == "medium" or line[headerIndex.index("confidence")] == "high":
            gene1_pos = line[headerIndex.index("breakpoint1")] # get exons how?
            gene2_pos = line[headerIndex.index("breakpoint2")]
            fusion_type = ""
            if line[headerIndex.index("#gene1")] == "BCR" and line[headerIndex.index("gene2")] == "ABL1":
                if int(gene1_pos.split(":")[1]) <= 23627388:
                    fusion_type = "Minor"
                elif int(gene1_pos.split(":")[1]) >= 23629346 and int(gene1_pos.split(":")[1]) <= 23637342 :
                    fusion_type = "Major"
                elif int(gene1_pos.split(":")[1]) >= 23651611:
                    fusion_type = " Micro"

            outline = [line[headerIndex.index("#gene1")], line[headerIndex.index("gene2")], gene1_pos, gene2_pos, fusion_type,
                    line[headerIndex.index("type")], line[headerIndex.index("split_reads1")],
                    line[headerIndex.index("split_reads2")], line[headerIndex.index("discordant_mates")],
                    line[headerIndex.index("coverage1")], line[headerIndex.index("coverage2")],
                    line[headerIndex.index("confidence")]]
            fusion_lines.append(outline)




## Xlsx file
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheetOver = workbook.add_worksheet("Overview")
worksheetBran = workbook.add_worksheet("Branford variants")
worksheetSNV = workbook.add_worksheet("SNV variants")
worksheetFusion = workbook.add_worksheet("Fusion")

headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
lineFormat = workbook.add_format({'top': 1})
tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})
textwrapFormat = workbook.add_format({'text_wrap': True})
italicFormat = workbook.add_format({'italic': True})
boldFormat = workbook.add_format({'bold': True})
redFormat = workbook.add_format({'font_color': 'red'})

''' Overview sheet '''
worksheetOver.write(0, 0, sample, headingFormat)
worksheetOver.write(1, 0, "Processing date: "+today.strftime("%B %d, %Y"))
worksheetOver.write_row(2, 0, emptyList, lineFormat)

worksheetOver.write(3, 0, "Created by: ")
worksheetOver.write(3, 4, "Valid from: ")
worksheetOver.write(4, 0, "Signed by: ")
worksheetOver.write(4, 4, "Document nr: ")
worksheetOver.write_row(5, 0, emptyList, lineFormat)

worksheetOver.write(6, 0, "Sheets:", tableHeadFormat)
worksheetOver.write_url(7, 0, "internal:'Branford variants'!A1", string="Branford variants")
worksheetOver.write_url(8, 0, "internal:'SNV variants'!A1", string="ABL1 variants")
worksheetOver.write_url(9, 0, "internal:'Fusions'!A1", string="Fusions")
worksheetOver.write_row(11, 0, emptyList, lineFormat)

worksheetOver.write(14, 0 , "Branford list used: "+str(snakemake.input.branford))
worksheetOver.write(15, 0, "Bedfile used for variantcalling: "+str(snakemake.input.bed) )


''' Branford list sheet '''
worksheetBran.write("A1", "Branford variants found", headingFormat)
worksheetBran.write("A3", "Sample: "+str(sample))
worksheetBran.write("A4", "Branford file: "+str(snakemake.input.branford))
worksheetBran.set_column('B:B', 28)
worksheetBran.set_column('D:D', 10)

worksheetBran.write_row("A6", branfordHeader, tableHeadFormat)
row = 6
col = 0
for line in branfordSNV:
    worksheetBran.write_row(row, col, line)
    row += 1


''' SNV ABL1 variants sheet '''
worksheetSNV.write("A1", "Variants found in ABL1 or BCR", headingFormat)
worksheetSNV.write("A3", "Sample: " + str(sample))
worksheetSNV.write("A4", "Bedfile used: " + str(snakemake.input.bed))
worksheetSNV.set_column('C:C', 10)

worksheetSNV.write_row("A6", allHeader, tableHeadFormat)
row = 6
for line in allSNVs:
    worksheetSNV.write_row(row, col, line)
    row += 1

''' Fusion sheet '''
worksheetFusion.set_column('C:D', 15)
worksheetFusion.set_column('F:F', 15)
worksheetFusion.write("A1", "Fusions detected with Arriba", headingFormat)
worksheetFusion.write("A3", "Sample: " + str(sample))
worksheetFusion.write("A5", "Primary fusion call from Arriba")

worksheetFusion.insert_image("A6", snakemake.input.jpg, {'x_scale': 0.25, 'y_scale': 0.25})
worksheetFusion.write_row("A40", fusionHeading, tableHeadFormat)
row = 40
for line in fusion_lines:
    worksheetFusion.write_row(row, col, line)
    row += 1

workbook.close()
