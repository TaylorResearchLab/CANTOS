import requests
import csv
import pandas as pd
import warnings

# Replace 'your_api_key' with your actual OpenAI API key
api_key = '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'

# Define the URL and headers for the API request
url = 'https://api.openai.com/v1/embeddings'
headers = {
    'Content-Type': 'application/json',
    'Authorization': f'Bearer {api_key}'
}

warnings.filterwarnings("ignore")

# Open the CSV file containing disease names
csv_file = 'diseases_8Nov23.csv'  # Replace with the path to your CSV file
encoding = 'utf-8'  # Try different encodings if needed

# Create an empty DataFrame to store embeddings
embeddings_df = pd.DataFrame(columns=['Disease'] + [f'embedding_{i}' for i in range(1500)])  # Initialize with columns

# Open the CSV file and read disease names from it
with open(csv_file, 'r', encoding=encoding) as file:
    csv_reader = csv.reader(file)
    #next(csv_reader)  # Skip the header row if it exists

    for row in csv_reader:
        disease = row[0]  # Assuming the disease name is in the first column
        print(disease)

        # Create the request data as a dictionary
        data = {
            'input': disease,
            'model': 'text-embedding-ada-002'
        }

        # Send a POST request to the OpenAI API
        response = requests.post(url, json=data, headers=headers, verify=False)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            result = response.json()
            embedding = result["data"][0]["embedding"]

            # Create a dictionary with disease name and embedding
            embedding_dict = {
                'Disease': disease,
                **{f'embedding_{i}': val for i, val in enumerate(embedding)}
            }

            # Add a new row to the DataFrame
            embeddings_df = embeddings_df.append(embedding_dict, ignore_index=True)
        else:
            print(f'Error for disease {disease}: Status Code {response.status_code}')
            print(f'Response: {response.text}')

# Save the DataFrame to a CSV file if needed
embeddings_df.to_csv('disease_embeddings.csv', index=False)

# Print the first few rows of the DataFrame
print(embeddings_df.head())
