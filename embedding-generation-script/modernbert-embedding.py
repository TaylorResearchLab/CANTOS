#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from transformers import AutoTokenizer, AutoModel
import torch
import numpy as np



# In[2]:


# Load ModernBERT model and tokenizer
model_name = "answerdotai/ModernBERT-large"  # Replace with the correct model name
tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
model.eval()  # Set model to evaluation mode



# In[3]:


# File paths
input_file = "biobert_embedding_terms.csv"  # Input CSV file
output_file = "mordernbert_embeddings.csv"  # Output CSV file

# Load tumor names from the CSV file
data = pd.read_csv(input_file)
if "Tumor_Names" not in data.columns:
    raise ValueError("The column 'Tumor_Names' is not present in the CSV file.")
tumor_names = data["Tumor_Names"].dropna().tolist()





# In[4]:


len(tumor_names)


# In[5]:


# Define batch size
batch_size = 5

def get_batch_embeddings(batch_texts, tokenizer, model):
    inputs = tokenizer(batch_texts, padding=True, truncation=True, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs)
    # Mean pooling to get a single vector per input
    attention_mask = inputs["attention_mask"]
    mask_expanded = attention_mask.unsqueeze(-1).expand(outputs.last_hidden_state.size())
    sum_embeddings = torch.sum(outputs.last_hidden_state * mask_expanded, 1)
    sum_mask = torch.clamp(mask_expanded.sum(1), min=1e-9)
    return sum_embeddings / sum_mask


# In[6]:


hidden_size = model.config.hidden_size
hidden_size


# In[7]:


# Initialize results DataFrame
columns = ["Tumor_Name"] + [f"Dim_{i}" for i in range(hidden_size)]  # Adjust dimension based on model
results_df = pd.DataFrame(columns=columns)


# In[9]:


# Process tumor names in batches
for i in range(0, len(tumor_names), batch_size):
    batch = tumor_names[i:i + batch_size]
    embeddings = get_batch_embeddings(batch, tokenizer, model)
    embeddings = embeddings.numpy()  # Convert to numpy array
    # Prepare batch results
    batch_results = pd.DataFrame(
        data=np.hstack([np.array(batch).reshape(-1, 1), embeddings]),
        columns=columns
    )
    # Append to results DataFrame
    results_df = pd.concat([results_df, batch_results], ignore_index=True)
    print(f"Processed batch {i // batch_size + 1}/{(len(tumor_names) + batch_size - 1) // batch_size}")


# In[12]:


# Convert results to DataFrame
results_df.shape



# In[13]:


# Save results incrementally to a CSV file
results_df.to_csv(output_file, index=False)
print(f"Embeddings saved to {output_file}")


# In[ ]:




