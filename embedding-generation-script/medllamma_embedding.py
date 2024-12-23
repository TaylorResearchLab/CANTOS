# Running on Respublica
import pandas as pd
import transformers
import torch
from transformers import AutoTokenizer, AutoModel

model_name = "johnsnowlabs/JSL-MedLlama-3-8B-v2.0"  # Replace with your MedLLaMA model path or Hugging Face model name
tokenizer = AutoTokenizer.from_pretrained(model_name)

model = AutoModel.from_pretrained(model_name)

csv_file = "biobert_embedding_terms.csv"  # Replace with your CSV file path
df = pd.read_csv(csv_file)

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

text_column = "Tumor_Names"  # Replace with your column name
texts = df[text_column].fillna("").tolist()  # Handle missing values as empty strings

embeddings = get_text_embeddings(texts)
df["embeddings"] = embeddings

# Save the DataFrame with embeddings to a new CSV file
output_file = "text_with_embeddings.csv"
df.to_csv(output_file, index=False)

