##############
#Written by Robin Sleith Fall 2019
#rsleith@smith.edu
#Version 1.1
#

#Step 0: Make a directory on Desktop named vizbin_centrifuge
#Step 1: Assemble your dataset and filter out contigs <500bp
#I use reformat.sh from bbtools in the terminal to do this
#reformat.sh in=assembly.fasta out=500_assembly.fasta minlength=500 -Xmx2g

#Step 2: Run centrifuge on this file (you could also use BLAST but it will take much longer...)
#centrifuge -p 20 -x nt -f 500_assembly.fasta -S 500_assembly.tab --report-file 500_assembly.tsv

#Step 3: Download minlength 500 assembly file to ~/Desktop/vizbin_centrifuge
#Get contig/scaffold IDs from assembly file. Copy the line below into textwrangler/bbedit, replace LKH with LKH number and copy the line to terminal, hit enter:
#cat ~/Desktop/vizbin_centrifuge/500_LKH_scaffolds.fasta |grep '>'|awk -F'>' '{print $2}'>~/Desktop/vizbin_centrifuge/left_ids_500_LKH_scaffolds.fasta

#Step 4: Download centrifuge file (ends in .tab) to ~/Desktop/vizbin_centrifuge
#Get top centrifuge result for each contig and clean up results (remove scaffolds classified as "species", "family" etc...)
#To do this copy the line below into textwrangler/bbedit and replace LKH with LKH number at beginning and end of the line!
#awk '!x[$1]++' ~/Desktop/vizbin_centrifuge/nt_LKH_scaffolds.fasta.tab|awk -F '\t' '$2!="species"' |awk -F '\t' '$2!="genus"'|awk -F '\t' '$2!="family"'|awk -F '\t' '$2!="order"'|awk -F '\t' '$2!="class"'|awk -F '\t' '$2!="phylum"'|awk -F '\t' '$2!="kingdom"'|awk -F '\t' '$2!="superkingdom"'|awk -F '\t' '$2!="no rank"' |awk -F '\t' '$2!="seqID"'  >~/Desktop/vizbin_centrifuge/clean_uniq_nt_LKH_scaffolds.fasta.tab

#Step 5: Download ncbi taxonomy table from github https://github.com/zyxue/ncbitax2lin and place in ~/Desktop/vizbin_centrifuge (ONLY DO THIS ONCE PER COMPUTER)
#read in the taxonomy table, ONLY DO THIS ONCE PER R SESSION
lineages <- read.csv("~/Desktop/vizbin_centrifuge/lineages-2019-02-20.csv", header=T, sep=',')
head(lineages)

#Step 6: Edit and run the following R script to add your files and vizbin output. I suggest running and debugging as you go
left <- read.table("~/Desktop/vizbin_centrifuge/left_ids_500_LKH_scaffolds.fasta")
left$id  <- 1:nrow(left)
head(left)
cent <- read.table("~/Desktop/vizbin_centrifuge/clean_uniq_nt_LKH_scaffolds.fasta.tab", header=F, sep='\t')
head(cent)
nrow(cent)
#below are two methods to subset centrifuge results, first by score, second by hit length relative to contig length
strong <- subset(cent, cent$score >50000)
strong <- subset(cent, cent$V6/cent$V7 >0.01)
nrow(strong)
head(strong)
#merge the centrifuge and taxonomy tables
#merged_lineages <- merge(x=strong, y=lineages, by.x="V3", by.y="tax_id", all=F) #uncomment this line to use subset of centrifuge results
merged_lineages <- merge(x=cent, y=lineages, by.x="V3", by.y="tax_id", all=F)
head(merged_lineages)
#merge the centrifuge/taxonomy with the contigs (neccessary as not every contig is classified by centrifuge)
merged_vizbin <- merge(x=left, y=merged_lineages,by.x = "V1", by.y = "V1" , all=T)
head(merged_vizbin)
#get in correct order for vizbin matching
final_merged_vizbin <- merged_vizbin[order(merged_vizbin$id), ]
head(final_merged_vizbin)

#Step 7: Run the vizbin app http://claczny.github.io/VizBin/ to get the path to the points file (About -> Show application log) 
#and copy to working directory (looks like /var/folders/ph/lzqf0q6s15s5n7bkhpdp3y200000gs/T/map6309370529122440509/points.txt) by typing the following in textwrangler/bbedit, changing LKH to LKH number and copying to terminal
#cp /var/folders/ph/lzqf0q6s15s5n7bkhpdp3y200000gs/T/map6309370529122440509/points.txt ~/Desktop/vizbin_centrifuge/500_LKH_points.txt
points <- read.table("~/Desktop/vizbin_centrifuge/500_LKH_points.txt", sep=",")
#plot(points)

#combine everything
tomap <- cbind(points, final_merged_vizbin)
head(tomap)
str(tomap)
##set up plotting
# Create new column filled with default colour
tomap$Color="blue"
# Set new column values to appropriate colours
tomap$Color[tomap$superkingdom=="Eukaryota"]="purple"
tomap$Color[tomap$superkingdom=="Bacteria"]="brown"
tomap$Color[tomap$superkingdom=="Archaea"]="pink"
tomap$Color[tomap$superkingdom=="Viruses"]="pink"
tomap$Color[tomap$kingdom=="Viridiplantae"]="green"
#tomap$Color[tomap$no.rank1=="Rhizaria"]="red"
tomap$Color[tomap$no.rank1=="Amoebozoa"]="red"
#tomap$Color[tomap$order =="Arcellinida"]="orange"
tomap$Color[tomap$kingdom=="Fungi"]="gold"
tomap$Color[tomap$kingdom=="Metazoa"]="black"

# Plot all points at once, using newly generated colours
plot(tomap$V2~tomap$V1, col=tomap$Color, pch=16, xlab="", ylab="")#xaxt='n',yaxt='n'
legend(17,-14, 
       legend = c("Archaea and Viruses", "Bacteria", "Viridiplantae", 
                  "Amoebozoa", "Fungi","Metazoa", "Other Eukaryota","Not_classified" ), 
       col = c("pink", "brown","green","red","gold","black", "purple","blue"),
       pch = c(16), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

#If you need to find a specific taxonomy value or level you will have to play around with the following lines 
#sort(unique(tomap$superkingdom))
#sort(unique(tomap$class))
#sort(unique(tomap$phylum))
#sort(unique(tomap$no.rank1))
#head(which(tomap == "Alveolata", arr.ind = TRUE))
#colnames(tomap[25])
#head(subset(tomap, tomap$no.rank1 =="Alveolata"))
##
#Below is an example for Rhizaria

