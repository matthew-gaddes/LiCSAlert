#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:10:36 2024

@author: matthew
"""

import pdb

#%%
def get_portal_public_volcanoes():
    """ Get the names of the volcanoes that are currently public on the 
    "COMET Volcanic and Magmatic Deformation Portal".  
    
    2024_06_11 | MEG | Written (with CGPT)
    """
    
    import requests
    from bs4 import BeautifulSoup
    from time import sleep
    
    # URL of the webpage to scrape
    url = "https://comet.nerc.ac.uk/comet-volcano-portal/volcano-index/Search-All"
    
    # Headers to mimic a real browser request (403 error otherwise)
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    
    # Use a session
    session = requests.Session()
    session.headers.update(headers)
    
    # Send a GET request to the webpage with retries
    for _ in range(3):  # Try up to 3 times
        response = session.get(url)
        if response.status_code == 200:
            break
        elif response.status_code == 403:
            print("403 Forbidden, retrying...")
            sleep(5)  # Wait for 5 seconds before retrying
        else:
            print(f"Failed to retrieve the webpage. Status code: {response.status_code}")
            break
    
    if response.status_code == 200:
        # Parse the webpage content
        soup = BeautifulSoup(response.content, "html.parser")
        
        # Find all the rows in the table (assuming each row represents a volcano)
        rows = soup.find_all("tr")
        
        # List to hold the names of the volcanoes
        volcano_names = []
    
        # Iterate through the rows and extract the volcano names
        for row in rows:
            # Find the column that contains the volcano name (assuming it's the first column)
            columns = row.find_all("td")
            if columns:
                volcano_name = columns[0].get_text(strip=True)
                volcano_names.append(volcano_name)
        
        # sort alphabetically
        volcano_names.sort()
            
        return volcano_names
    else:
        print(f"Failed to retrieve the webpage after retries. Status code: {response.status_code}")
        
        
#%%