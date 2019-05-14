#!/usr/bin/env python3
# Copyright 2018 Lael D. Barlow
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
"""Script for adding adding filenames as labels to the bottom left corner of
all pdf files in a given directory (outputs copies of the pdf files in a
subdirectory). This is useful for adding names to trees output by FigTree.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))

from PyPDF2 import PdfFileWriter, PdfFileReader
import io
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

import glob

cmdln = sys.argv
dirpath = cmdln[1]

# Get paths for pdfs to be labelled.
filenames = glob.glob(os.path.join(dirpath, '*.pdf'))
outdir = os.path.join(dirpath, 'labelled_pdfs')
os.mkdir(outdir)

for filename in filenames:
    # Define output filename.
    outfilename = os.path.join(outdir, os.path.basename(filename))

    packet = io.BytesIO()
    # create a new PDF with Reportlab
    can = canvas.Canvas(packet, pagesize=letter)
    can.drawString(30, 30, os.path.basename(filename))
    can.save()

    #move to the beginning of the StringIO buffer
    packet.seek(0)
    new_pdf = PdfFileReader(packet)
    # read your existing PDF
    existing_pdf = PdfFileReader(open(filename, "rb"))
    output = PdfFileWriter()
    # add the "watermark" (which is the new pdf) on the existing page
    page = existing_pdf.getPage(0)
    page.mergePage(new_pdf.getPage(0))
    output.addPage(page)
    # finally, write "output" to a real file
    outputStream = open(outfilename, "wb")
    output.write(outputStream)
    outputStream.close()
