import pandas as pd

if __name__ == "__main__":
    filename = input("enter filename\n")
    df = pd.read_csv(filename)
    print("Shape: ")
    print(df.shape)

    print("Columns: ")
    print(df.columns)

    print("First 5 rows: ")
    print(df.head())

    print("Last 5 rows: ")
    print(df.tail())