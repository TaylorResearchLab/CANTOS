#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install sentencepiece')
get_ipython().system('pip install -U accelerate')


# In[ ]:





# In[2]:


import pandas as pd
import torch
from transformers import LlamaTokenizer, LlamaModel
import sentencepiece


# In[3]:


MODEL_NAME = "chaoyi-wu/MedLLaMA_13B"


# In[4]:


tokenizer = LlamaTokenizer.from_pretrained(MODEL_NAME)
model = LlamaModel.from_pretrained(MODEL_NAME, device_map="auto")


# In[5]:


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)  # Move the model to the selected device


# In[6]:


# Load the CSV file
csv_file = "biobert_embedding_terms.csv"  # Replace with your CSV file path
data = pd.read_csv(csv_file)


# In[7]:


data.head(10)


# In[8]:


def generate_embeddings(text):
    inputs = tokenizer(text, return_tensors="pt", truncation=True).to(device)  # Move inputs to the correct device
    with torch.no_grad():
        outputs = model(**inputs)
    # Use the hidden states of the last layer as embeddings
    return outputs.last_hidden_state.mean(dim=1).squeeze().cpu().numpy()



# In[10]:


embeddings = []
iter=1
for word in data['Tumor_Names']:
    embedding = generate_embeddings(word)
    embeddings.append(embedding)
    if iter % 100 == 0:
        print(iter)
    iter=iter+1  


# In[11]:


embedding_dim = embeddings[0].shape[0]
embedding_columns = [f"embedding_{i}" for i in range(embedding_dim)]
embedding_df = pd.DataFrame(embeddings, columns=embedding_columns)
result_df = pd.concat([data, embedding_df], axis=1)


# In[12]:


result_df.head(10)


# In[13]:


output_csv_file = "medllama-13b.csv"
result_df.to_csv(output_csv_file, index=False)


# In[ ]:


iter=1
for Tumor_Names in data['Tumor_Names']:
    embedding = generate_embeddings(Tumor_Names)
    embeddings.append(embedding)
    if iter % 100 == 0:
        print(iter)
    iter=iter+1        


# In[ ]:


def get_text_embeddings(texts):
    embeddings = []
    iter=1
    for text in texts:
        inputs = tokenizer(text, return_tensors="pt", padding=True, truncation=True)
        with torch.no_grad():
            outputs = model(**inputs)
        # Compute the average embedding across tokens (sequence level embedding)
        embedding = outputs.last_hidden_state.mean(dim=1).squeeze().numpy()
        embeddings.append(embedding)
        iter=iter+1
        if iter % 100 == 0:
            print(iter)
    return embeddings


# In[ ]:


text_column = "Tumor_Names"  # Replace with your column name
texts = data[text_column].fillna("").tolist()  # Handle missing values as empty strings


# In[ ]:


embeddings = get_text_embeddings(texts)


# In[ ]:




