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
"""For sending emails from scripts. Some code taken from
http://naelshiab.com/tutorial-send-email-python/

This works for sending emails from gmail accounts for which "less secure apps"
have been enabled in the security settings.  
"""
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import settings

def send_eml(subject, body):
    """Sends an email using the from and to addresses and password defined in
    the settings.py file.
    """
    # Compose email.
    fromaddr = settings.fromaddr
    pswd = settings.pswd
    toaddr = settings.toaddr

    # Send email.
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = subject
    msg.attach(MIMEText(body, 'plain'))
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(fromaddr, pswd)
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()

if __name__ == '__main__':
    send_eml('Test subject', 'Test body')
