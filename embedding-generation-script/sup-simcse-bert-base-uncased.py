#!/usr/bin/env python
# coding: utf-8

# In[1]:


from transformers import AutoTokenizer, AutoModel
import torch
import pandas as pd


# In[2]:


# Load the pretrained SimCSE model and tokenizer
model_name = "princeton-nlp/sup-simcse-bert-base-uncased"  # Supervised SimCSE model
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)



# In[3]:


# Function to generate embeddings for a single text
def generate_embedding(text):
    # Tokenize the text
    inputs = tokenizer(text, return_tensors="pt", truncation=True, padding=True, max_length=512)
    # Forward pass through the model
    with torch.no_grad():
        outputs = model(**inputs)
    # Extract the [CLS] token embedding
    cls_embedding = outputs.last_hidden_state[:, 0, :]
    return cls_embedding.squeeze().numpy()


# In[4]:


# Load the CSV file
csv_file = "biobert_embedding_terms.csv"  # Replace with your file path
data = pd.read_csv(csv_file)

# Specify the column containing text
text_column = "Tumor_Names"  # Replace with the name of your column

# Check if the column exists
if text_column not in data.columns:
    raise ValueError(f"Column '{text_column}' not found in the CSV file.")



# In[5]:


# Generate embeddings for each row
embeddings = []
iter=1
for text in data[text_column]:
    embedding = generate_embedding(text)
    embeddings.append(embedding)
    if iter%1000==0:
        print(iter)
    iter=iter+1    


# In[6]:


embedding_dim = embeddings[0].shape[0]
embedding_columns = [f"embedding_{i}" for i in range(embedding_dim)]
embedding_df = pd.DataFrame(embeddings, columns=embedding_columns)
result_df = pd.concat([data, embedding_df], axis=1)


# In[7]:


result_df.head(10)


# In[8]:


# Save embeddings to a new CSV file
output_file = "sup-simcse-bert-base-uncased-embeddings.csv"
result_df.to_csv(output_file, index=False)
print(f"Embeddings saved to {output_file}")


# In[ ]:




