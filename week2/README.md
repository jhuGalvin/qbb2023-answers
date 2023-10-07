## Exercise 1: Coverage Simulator

## Step 1.1

(1Mb x 3)/100bp = 30,000 reads for 3x coverage

## Step 1.3

In my simulation at 3x coverage, 5.18% of the genome does not have any coverage. This matches the expected coverage under the Poisson distribution at lambda = 3 as well as the normal distribution, although the data is skewed slightly to the left of the normal distribution. 

## Step 1.4

In my simulation at 10x coverage, 0.0065% of the genome does not have any coverage. This matches the expected coverage under the Poisson distribution at lambda = 10 as well as the normal distribution, although the data is skewed slightly to the left of the normal distribution. 

## Step 1.5

In my simulation at 30x coverage, 0.0001% of the genome does not have any coverage. This matches the expected coverage under the Poisson distribution at lambda = 30 as well as the normal distribution, although the data is skewed slightly to the left of the normal distribution. 


## Exercise 2: de Bruijn graph construction
## Step 2.4

Text used in terminal to create graph
(base) [~/qbb2023-answers/week2 $]conda activate graphviz
(graphviz) [~/qbb2023-answers/week2 $]dot -Tpng edge.txt -o ex2_digraph.png

## Step 2.5

ATTCATTGATTCTTATTCTTAT

## Step 2.6
To reconstruct the genome, ideally I would have a reference to map my reads to and deep coverage throughout my libraries. Regions of structural variation or high copy number/repeats would be more difficult to construct, so access to long read sequencing would also be nice.

## Step 3.1

Prompt: Hello friend; please write an essay telling me why a young scientist or student might want to learn about genomics. Please include references as well. Thanks!


Title: The Genomic Frontier: A Journey for Young Scientists and Students

Introduction

Genomics, the study of an organism's complete set of genes, DNA sequences, and their interactions, represents a captivating and dynamic field in the realm of biology and genetics. For young scientists and students, the pursuit of genomics education offers a gateway to a world filled with innovation, scientific discovery, and a host of practical applications that have far-reaching implications. In this essay, we will explore the compelling reasons why a young scientist or student might want to learn about genomics, supported by references and examples that underline the significance of this field.

Advancing Scientific Knowledge
Genomics is at the forefront of modern biology and genetics research. It seeks to unravel the fundamental mysteries of life by decoding the genetic blueprint of organisms. A foundational understanding of genomics empowers young scientists to contribute to cutting-edge research and discovery. By studying the genomic basis of life, one can delve into the intricacies of evolution, species diversity, and the molecular mechanisms underpinning various biological processes (Venter et al., 2001).

Medical Applications
Genomics is a cornerstone of the healthcare industry. For students aspiring to become healthcare professionals, the study of genomics is indispensable. Genomic medicine enables the development of personalized treatment plans based on an individual's genetic makeup, promising more effective and targeted therapies. It has been pivotal in diagnosing genetic disorders, tailoring cancer treatments, and understanding the hereditary basis of diseases (Collins et al., 2003).

Agricultural Advancements
In an era marked by climate change and a growing global population, genomics is revolutionizing agriculture. The knowledge of genomics can be harnessed to develop more resilient and productive crops, enhance livestock breeds, and combat agricultural challenges such as pests and diseases. Students passionate about agriculture can drive innovation and contribute to food security and sustainability (Varshney et al., 2015).

Biotechnology and Biopharmaceuticals
The biotechnology and biopharmaceutical industries rely heavily on genomics. The production of genetically modified organisms (GMOs), gene therapies, and the design of therapeutic drugs are just a few examples of how genomics is transforming these sectors. Young scientists who study genomics can spearhead research into groundbreaking therapies, vaccines, and bioengineering solutions (Lander et al., 2001).

Conservation Biology
For the environmentally conscious and conservation-minded, genomics plays a vital role in preserving biodiversity and managing endangered species. It can provide insights into genetic diversity, aid in combating illegal wildlife trafficking, and support habitat preservation efforts. Genomics equips students with the tools to protect and conserve the planet's rich ecosystems (Morin et al., 2010).

Career Opportunities
A profound understanding of genomics opens up a wide array of career opportunities. Graduates with genomics expertise can pursue careers in academic research, healthcare, biotechnology, forensics, or bioinformatics, among others. The field's diversity and the growing demand for genomics specialists make it an attractive avenue for young scientists (Lewin et al., 2004).

Ethical and Societal Implications
Studying genomics also requires grappling with the ethical, legal, and societal implications of genetic research and technology. This knowledge is critical for shaping responsible policies, regulations, and addressing the ethical dilemmas associated with genetic advancements. Students and scientists who engage with these ethical questions contribute to responsible innovation (Knoppers et al., 2016).

Conclusion

In conclusion, genomics is a captivating field with an extensive range of reasons why a young scientist or student should consider learning about it. From advancing scientific knowledge to enabling medical breakthroughs, from driving agricultural innovations to propelling careers in diverse industries, genomics holds the promise of making meaningful contributions to science, medicine, and society as a whole. It is a frontier of knowledge that beckons the curious and dedicated minds of the next generation.

References:

Venter, J. C., et al. (2001). The Sequence of the Human Genome. Science, 291(5507), 1304-1351.
Collins, F. S., et al. (2003). A Vision for the Future of Genomic Research. Nature, 422(6934), 835-847.
Varshney, R. K., et al. (2015). Genomics‐assisted breeding for crop improvement. Trends in Plant Science, 20(12), 621-631.
Lander, E. S., et al. (2001). Initial sequencing and analysis of the human genome. Nature, 409(6822), 860-921.
Morin, P. A., et al. (2010). Application of SNP genotyping in conservation genetics and management of endangered species. Molecular Ecology Resources, 10(1), 83-92.
Lewin, H. A., et al. (2004). The first international conference on livestock genomics. Genome Research, 14(2), 233-238.
Knoppers, B. M., et al. (2016). Legal uncertainties arising from direct-to-consumer personal genome testing. Journal of Law and the Biosciences, 3(3), 619-662.

## Step 3.2

ChatGPT wrote an essay that makes logical sense but lacks flow and uses a variety of generalizing phrases. The essay seems more like a list one might see on a website than a true essay, but is useful for generating ideas. The citations are real for the most part, some seem difficult to track down or may not really exist. The sources do not always seem to actually have content ChatGPT is using it for, like the Science articles about the sequence of the human genome.