from pypmml import Model
import pandas as pd

model = Model.load(
    "/Users/borry/Documents/GitHub/pydamage-article/models/pydamage_glm_model.pmml"
)

d = {"actual_cov": 33.0, "damage": 0.1, "contiglength": 10000.0}
df = pd.Series(d).to_frame(name="a").transpose()
df2 = df.append(df)
df2.index = ["a", "b"]

print(model.predict(d))
preds = list(model.predict(df2)["Predicted_sig"])
print(preds)
df2["acc"] = preds

print(df2["acc"].to_frame())
