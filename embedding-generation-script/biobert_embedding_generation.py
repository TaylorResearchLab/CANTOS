import pandas as pd
from transformers import AutoTokenizer, AutoModel
import torch

# Step 1: Read the CSV file into a DataFrame
csv_file = "biobert_embedding_terms.csv"  # Replace with your file path
df = pd.read_csv(csv_file)

# Step 2: Select the text column for generating embeddings
# Replace 'text_column' with the actual column name containing text
text_data = df['Tumor_Names'].dropna().tolist()

# Step 3: Load the BioBERT tokenizer and model
model_name = "dmis-lab/biobert-base-cased-v1.1"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)

# Step 4: Function to generate embeddings and format them into a DataFrame
def generate_embeddings(texts):
    all_embeddings = []
    iter=1
    for text in texts:
        # Tokenize the input text
        inputs = tokenizer(text, return_tensors="pt", truncation=True, padding=True, max_length=512)
        if iter%1000==0:
          print(iter)
        iter=iter+1
        # Get the model's output
        with torch.no_grad():
            outputs = model(**inputs)
        
        # Extract the embedding (using the [CLS] token representation)
        cls_embedding = outputs.last_hidden_state[:, 0, :].squeeze(0).numpy()
        all_embeddings.append(cls_embedding)
    
    # Convert the list of embeddings into a DataFrame (each dimension is a column)
    embedding_df = pd.DataFrame(all_embeddings, columns=[f"dim_{i+1}" for i in range(cls_embedding.shape[0])])
    return embedding_df


# Step 5: Generate embeddings for the text data

embeddings_df = generate_embeddings(text_data)

# Step 6: Combine the original DataFrame with the embeddings
# Resetting index to ensure alignment if the text column had NaNs removed
df = df.reset_index(drop=True)
final_df = pd.concat([df, embeddings_df], axis=1)

# Save the final DataFrame with embeddings to a CSV file
final_df.to_csv("output_with_embeddings.csv", index=False)