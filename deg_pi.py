#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# To run: streamlit run <<program name>>>
# Then go to http://localhost:8501/ with browser

# Date: 26 October 2022

""" Streamlit demo to open an sqlite3 database """

import streamlit as st
import pandas as pd
import plotly.express as px
from PIL import Image

# DB Mgmt
import csv
import re
import sqlite3
import subprocess
import os
import math
import dash_bio

def file_selector(folder_path='.'): ## pick and load a database from the local directory
    filenames = os.listdir(folder_path)
    selected_filename = st.selectbox('Select a file', filenames)
    return os.path.join(folder_path, selected_filename)

# Fxn Make Execution
def sql_executor(myCommand_str):
    """ function to complete the query and parse results from a query."""
    c.execute(myCommand_str)
    data = c.fetchall()
    return data
# end of sql_executor()

# check and open file
def get_file():
    
    columns = ["B","t"]
    data = pd.read_csv("./r.csv")
    data

    df = pd.read_csv("./r.csv", usecols=columns)
    x = df.t
    y = df.B8

    return x,y

# writes the keyword into a file to be read by r
def add_data_to_file(file_name,word):
    file_object = open(file_name, "w")
    file_object.write(word + "\n")
    file_object.close()

# def get_genes_and_score():
#     with open("./r.csv", 'r') as file:
#         csvreader = csv.reader(file)
#         is_header = True
#         samples_not_included = 0 #because gene symbol is unavailable 


#         with open('sig_genes.csv', 'w', newline='') as file:
#             writer = csv.writer(file)
            
#             for row in csvreader:
#                 # does not try to see if header is significant
#                 if is_header:
#                     is_header = False
#                     writer.writerow(row)
#                 else: 
#                     p_val = float(row[1])
#                     logFC = float(row[5])
#                     gene = row[6]

#                     if gene != "":  # a gene symbol is available
#                         if p_val <= 0.05 and logFC < 0.5: #difference is statiscally significant
#                             writer.writerow(row)
#                     else:
#                         samples_not_included += 1 # excluded because I don't have the gene symbol or name

#         genes = ""
#         score = ""
#         is_header2 = True
#         with open('sig_genes.csv', 'r', newline='') as file:   
#             for row in csv.reader(file):
#                 if is_header2:
#                     is_header2 = False
#                 else:
#                     genes += row[6] + ","
#                     score += str(math.log(float(row[1]),10)) + ","
#         genes = genes[:-1]
#         score = score[:-1]

#         add_data_to_file("gene_list_Webgestalt.txt", genes)
#         add_data_to_file("score_list_Webgstalt.txt", score)


def display_webgestalt_results():	
    image = Image.open('./Project_test/goslim_summary_test.png')

    st.image(image, caption=None, width=None, use_column_width=None, clamp=False, channels="RGB", output_format="auto")
    
    # image 1 
    # Size of the image in pixels (size of original image)
    # (This is not mandatory)
    width, height = image.size
    
    # Setting the points for cropped image
    left = 0
    top = 0
    right = width / 3 
    bottom = height
    
    # Cropped image of above dimension
    im1 = image.crop((left, top, right, bottom))

    # image 2
    # Size of the image in pixels (size of original image)
    # (This is not mandatory)
    width, height = image.size
    
    # Setting the points for cropped image
    left = width / 3
    top = 0
    right = 2 * width / 3 
    bottom = height
    
    # Cropped image of above dimension
    im2 = image.crop((left, top, right, bottom))

    # image 3
    # Size of the image in pixels (size of original image)
    # (This is not mandatory)
    width, height = image.size
    
    # Setting the points for cropped image
    left = 2 * width / 3
    top = 0
    right = width
    bottom = height
    
    # Cropped image of above dimension
    im3 = image.crop((left, top, right, bottom))

    st.image(im1, caption=None, width= None, use_column_width=None, clamp=False, channels="RGB", output_format="auto")
    st.image(im2, caption=None, width= None, use_column_width=None, clamp=False, channels="RGB", output_format="auto")
    st.image(im3, caption=None, width= None, use_column_width=None, clamp=False, channels="RGB", output_format="auto")

    pathway_overview = pd.read_csv('./Project_test/enrichment_results_test.csv',usecols=[0,1,4,5])
    st.write(pathway_overview)

    with open("./Project_test/enrichment_results_test.txt", 'r') as myfile: 
        with open("./Project_test/enrichment_results_test.csv", 'w') as csv_file:
            for line in myfile:

                # Replace every tab with comma
                fileContent = re.sub("\t", ",", line)
       
                # Writing into csv file
                csv_file.write(fileContent)
        
        myfile.close()
        csv_file.close()


        results = open("./Project_test/enrichment_results_test.csv", 'r')
        result_reader = csv.reader(results)
        st.title("Enriched Pathways")



        for line in result_reader:
            if line[0] != "geneSet":
                
                st.header(line[0] + ": " + line[1])
                st.subheader("Enrichment Score: " + str("%.4f" % (float(line[3]))))
                st.subheader("P-value: " + str("%.4f" % (float(line[5]))))

                imagetitle = './Project_test/Project_test_GSEA/' + line[0] + '.png'
                image = Image.open(imagetitle)
                st.image(image, caption=None, width= None, use_column_width=None, clamp=False, channels="RGB", output_format="auto")


                st.subheader("Genes Involved:") 


                gene_string = line[11].replace(";", ", ")
                
                st.caption(gene_string) #print out string of genes in pathway

                st.header("")


def create_volcano_plot():
    df = pd.read_csv('r.csv')

    fig = dash_bio.VolcanoPlot(
        dataframe = df,
        effect_size= 'logFC',
        p='adj.P.Val',
        snp='Gene.symbol',
        gene='Gene.title',
        annotation=None,
        logp=True,
        xlabel=None,
        ylabel='-log10(p)',
        point_size=2,
        col=None,
        effect_size_line=[-0.5,0.5],
        effect_size_line_color='yellow',
        effect_size_line_width=0.5,
        genomewideline_value = 1.30103,
        genomewideline_color='blue',
        genomewideline_width=0.1,
        highlight=True,
        highlight_color='red',
    )


    st.plotly_chart(fig, use_container_width=True)


def sort_into_groups():

    # categorize 

    # creates a new coloumn for the group number to go in the same file
    # creates a header to match the rest of the data
    csv_input = pd.read_csv('sample_type.csv')
    csv_input['group_number'] = ''


    # opens the file and determines line by line group type by looking for keywords in sample
    with open('sample_type.csv','r') as csvinput:
        row_num = 0
        result_string = "A"
        for row in csv.reader(csvinput):
            title = row[0].lower()
            if title != "title":
                if 'cancer' in title or 'tumor' in title: 
                    group = '1' #affect tissue
                elif 'surround' in title or 'control' in title or 'normal' in title:
                    group = '0' #unaffected tissue
                else: 
                    group = '1' #affected tissue    
                    
                result_string += group
            
                csv_input.loc[row_num,'group_number'] = group
            row_num += 1

        csv_input.to_csv('output.csv',index=False)

        file_object = open("group_membership.txt", "w")
        file_object.write(result_string + "\n")
        file_object.close()

def find_degs():

    with open("./r.csv", 'r') as file:
        csvreader = csv.reader(file)
        is_header = True
        samples_not_included = 0 #because gene symbol is unavailable

        with open('sig_genes.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            
            for row in csvreader:
                # does not try to see if header is significant
                if is_header: 
                    is_header = False
                else: 
                    p_val_adj = float(row[1])
                    p_val = float(row[2])
                    logFC = float(row[5])
                    gene = row[6]

                    if gene != "":  # a gene symbol is available
                        if p_val_adj <= 0.05 and logFC < 0.5 and logFC > -0.5: #difference is statiscally significant
                            writer.writerow(row)
                    else:
                        samples_not_included += 1 # excluded because I don't have the gene symbol or name

def find_pathways():
    genes = ""
    score = ""
    file_object = open("webgestaltData.txt", "w")
    with open('sig_genes.csv', 'r', newline='') as file:
        rowcount = 0
        for row in csv.reader(file):
            text = row[6] + "\t" + str("%.4f" % (float(row[5])))
            file_object.write(text + "\n")
    file_object.close() #close file
        
    txt_file = 'webgestaltData.txt'
    base = os.path.splitext(txt_file)[0]
    os.rename(txt_file, base + '.rnk')


def deg_output():
    st.header("Differentally Expressed Genes")

    deg_file = pd.read_csv("r.csv")
    st.write(deg_file)
    create_volcano_plot()

def display_about():
    st.header("What do the results mean?")
    with st.container():
            st.subheader("Variables")
            st.markdown("p_value: the raw probablity that the difference occured simply due to change\n\nadj_p_val: p-value after adjustment for multiple testing or false positives, generally accepted to be more accurate\n\nlogFC: Log2-fold change, fold change determines whether a gene is upregulated or downregulated and the respective magnitude")


    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Volcano Plot")
        st.markdown("A volcanot plot is a type of scatter plot that analyzes statistical significance (x-axis) versus magnitude of change (y-axis). The upregulated genes have a positive magnitude of change, and the more upregulated are farther away from 0. Similarly, the downregulated genes have a negative magnitude of change, and the more downregulated are farther away from 0. The most statistically signficant genes are towards the top of the graph.\n\nThe red color signifies the genes that are both statistically significant (adjusted p-value < 0.5) and have a fold value between -0.5 and 0.5. These are the DEGs. The blue represents the genes the do not meet the criteria and thus are not considered a DEG.")
    with col2: 
        st.text("\n\n\n\n")
        st.image("volcanoPlot.png")
    
def demo():
    st.title("DEG-PI")

    menu = ["Home","DEGs","Enriched Pathways","About"]
    choice = st.sidebar.selectbox("Menu",menu)

    if choice == "Home":
        st.header("Welcome to DEG-PI (Differentally Expressed Gene and Pathway Identifier):\n")
        st.subheader("A Computational Tool the Identifies Differentially Expressed Genes and Enriched Pathways in Cancer Tissue")
        with st.form(key='query_form'):
            keyword = st.text_area("To Get Started Enter a Keyword or GEO Accession")
            submit_button = st.form_submit_button("Submit")

        if submit_button:
            # check validity
            st.success("We are collecting your data. This could take up to 15 minutes.") # print up a message when you click "Execute"
            
                 # get keyword from user and add to file to be read by R program
            add_data_to_file("keyword.txt", keyword)
            process1 = subprocess.Popen(["Rscript", "getData.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            result1 = process1.communicate()
            
            sort_into_groups()

            # Use categories to get genes that are overly expressed
            process2 = subprocess.Popen(["Rscript", "/Users/abbywaryanka/seniorComp/Streamlit2/findsGenes.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            result2 = process2.communicate()

            # use data to find differentially expressed genes
            find_degs()

            #set up data for Webgetstalt
            find_pathways()

            #Run Webgestaltr
            process3 = subprocess.Popen(["Rscript", "find_common_pathways.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            result3 = process3.communicate()

            st.caption("Data has successfully been collected. Please navigate to respective pages to view results.")

    elif choice == "DEGs":

        deg_output()

                
    elif choice == "Enriched Pathways":
        st.header("Pathway Enrichment Anlaysis")
        display_webgestalt_results()

    elif choice == "About":
        display_about()


if __name__ == '__main__':
    demo()
