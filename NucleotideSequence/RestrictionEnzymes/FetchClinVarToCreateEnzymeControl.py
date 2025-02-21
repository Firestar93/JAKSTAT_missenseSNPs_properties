import requests
import time
import pandas as pd

full_output_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\ClinVarBenignMutations\\clinvar_significance.tsv"
benign_output_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\ClinVarBenignMutations\\clinvar_benign_variants.tsv"

extracted_variants_list = [
    ("STAT1", "His58Tyr", "rs751586208"),
    ("STAT1", "Arg241Gln", "rs146273341"),
    ("STAT1", "Ile265Val", "rs148775168"),
    ("STAT1", "Val266Ile", "rs41473544"),
    ("STAT1", "Val455Ile", "rs371982540"),
    ("STAT1", "Asn515Ser", "rs190269533"),
    ("STAT1", "Ala531Thr", "rs148573907"),
    ("STAT1", "Ala589Ser", "rs745491762"),
    ("STAT1", "Arg619Gln", "rs369060692"),
    ("STAT1", "Pro696His", "rs138723664"),
    ("STAT3", "Val461Leu", "rs149214040"),
    ("STAT3", "Ile498Val", "rs146620441"),
    ("STAT3", "Val507Phe", "rs145786768"),
    ("STAT3", "Gly618Arg", "rs2081548277"),
    ("STAT3", "Gly618Arg", "rs2081548277"),
    ("STAT3", "Gly618Arg", "rs2081548277"),
    ("STAT3", "Tyr614Phe", "rs769031989"),
    ("STAT3", "Asn647Ile", "rs770986654"),
    ("STAT3", "Asn647Ile", "rs770986654"),
    ("STAT3", "Asn647Ile", "rs770986654"),
    ("STAT3", "Asn647Ile", "rs770986654"),
    ("STAT3", "Asn647Ile", "rs770986654"),
    ("STAT3", "Asp661Tyr", "rs747639500"),
    ("STAT3", "Asp661Tyr", "rs747639500"),
    ("STAT3", "Asp661Tyr", "rs747639500"),
    ("STAT3", "Asp661Tyr", "rs747639500"),
    ("STAT3", "Asp661Tyr", "rs747639500"),
    ("STAT3", "Ala702Thr", "rs747667389"),
    ("STAT3", "Gly743Val", "rs151033214"),
    ("STAT3", "Ser763Leu", "rs140604473"),
    ("STAT4", "Ile115Val", "rs3024839"),
    ("STAT4", "Glu128Val", "rs140675301"),
    ("STAT4", "Val143Met", "rs370819441"),
    ("STAT4", "Arg240Gln", "rs61756200"),
    ("STAT4", "Asp256Asn", "rs771061948"),
    ("STAT4", "Leu269Ile", "rs35279173"),
    ("STAT4", "Thr298Ile", "rs200982266"),
    ("STAT4", "Ala361Thr", "rs751746803"),
    ("STAT4", "Gly393Ala", "rs199633613"),
    ("STAT4", "Thr446Ile", "rs141331848"),
    ("STAT4", "Asn475Ser", "rs146386562"),
    ("STAT4", "Asp513Val", "rs147007636"),
    ("STAT5B", "Leu94Phe", "rs199645527"),
    ("STAT5B", "Arg100Cys", "rs199894785"),
    ("STAT5B", "Ala130Val", "rs2277619"),
    ("STAT5B", "Glu315Ala", "rs572536541"),
    ("STAT5B", "Arg353His", "rs143171571"),
    ("STAT5B", "Arg353Cys", "rs762833594"),
    ("STAT5B", "Ala510Val", "rs200200711"),
    ("STAT5B", "Asn642His", "rs938448224"),
    ("STAT5B", "Asn642His", "rs938448224"),
    ("STAT5B", "Asn642His", "rs938448224"),
    ("STAT5B", "Asn642His", "rs938448224"),
    ("STAT5B", "Asn642His", "rs938448224"),
    ("JAK2", "Gly13Val", "rs759031245"),
    ("JAK2", "Ile19Val", "rs150159583"),
    ("JAK2", "Ser46Tyr", "rs138655335"),
    ("JAK2", "Gly48Glu", "rs143227399"),
    ("JAK2", "Leu113Val", "rs143103233"),
    ("JAK2", "Gly127Asp", "rs56118985"),
    ("JAK2", "Gln175Glu", "rs756128160"),
    ("JAK2", "Ile237Phe", "rs376125987"),
    ("JAK2", "Lys244Arg", "rs62637625"),
    ("JAK2", "Asn337Asp", "rs149683525"),
    ("JAK2", "Ile354Thr", "rs371907546"),
    ("JAK2", "Leu383Val", "rs143124074"),
    ("JAK2", "Val392Met", "rs200018153"),
    ("JAK2", "Leu393Val", "rs2230723"),
    ("JAK2", "Gly417Ser", "rs190968273"),
    ("JAK2", "Cys480Phe", "rs62637623"),
    ("JAK2", "Arg487Cys", "rs764423560"),
    ("JAK2", "Arg564Leu", "rs368927897"),
    ("JAK2", "Arg564Leu", "rs368927897"),
    ("JAK2", "Val567Ala", "rs587778408"),
    ("JAK2", "Gly571Ser", "rs139504737"),
    ("JAK2", "Gly571Ser", "rs139504737"),
    ("JAK2", "Gly571Ser", "rs139504737"),
    ("JAK2", "Gly571Ser", "rs139504737"),
    ("JAK2", "His587Asn", "rs149705816"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Phe", "rs77375493"),
    ("JAK2", "Val617Ile", "rs77375493"),
    ("JAK2", "Val617Ile", "rs77375493"),
    ("JAK2", "Val617Ile", "rs77375493"),
    ("JAK2", "Arg683Gly", "rs1057519721"),
    ("JAK2", "Arg683Gly", "rs1057519721"),
    ("JAK2", "Arg683Gly", "rs1057519721"),
    ("JAK2", "Arg683Gly", "rs1057519721"),
    ("JAK2", "Arg683Gly", "rs1057519721"),
    ("JAK2", "Arg683Gly", "rs1057519721"),
    ("JAK2", "Asn691His", "rs151160183"),
    ("JAK2", "Ile724Thr", "rs372254348"),
    ("JAK2", "Ser759Tyr", "rs766524586"),
    ("JAK2", "Ser797Cys", "rs201992086"),
    ("JAK2", "Arg839Gln", "rs747381013"),
    ("JAK2", "Ile899Thr", "rs200282557"),
    ("JAK2", "Glu890Gly", "rs368599778"),
    ("JAK2", "Asp894Gly", "rs757780497"),
    ("JAK2", "Leu892Val", "rs201551707"),
    ("JAK2", "Ile1051Thr", "rs375671491"),
    ("JAK2", "Lys1053Arg", "rs142094756"),
    ("JAK2", "Arg1063His", "rs41316003"),
    ("JAK2", "Arg1063His", "rs41316003"),
    ("JAK2", "Arg1063His", "rs41316003"),
    ("JAK2", "Arg1063His", "rs41316003"),
    ("JAK2", "Asn1108Ser", "rs142269166"),
    ("JAK2", "Asn1108Ser", "rs142269166"),
    ("JAK2", "Asn1108Ser", "rs142269166"),
    ("JAK3", "Thr8Met", "rs145500023"),
    ("JAK3", "Pro12Leu", "rs56061056"),
    ("JAK3", "Arg40His", "rs56384680"),
    ("JAK3", "Ile63Val", "rs144405201"),
    ("JAK3", "Val90Met", "rs1016346013"),
    ("JAK3", "Arg121His", "rs143586866"),
    ("JAK3", "Pro132Thr", "rs3212723"),
    ("JAK3", "Pro132Thr", "rs3212723"),
    ("JAK3", "Pro132Thr", "rs3212723"),
    ("JAK3", "Pro132Thr", "rs3212723"),
    ("JAK3", "Val217Met", "rs202167678"),
    ("JAK3", "Arg222His", "rs199868795"),
    ("JAK3", "Ala261Val", "rs777830679"),
    ("JAK3", "Val299Ala", "rs571404212"),
    ("JAK3", "Thr381Asn", "rs373046546"),
    ("JAK3", "Arg395His", "rs143038064"),
    ("JAK3", "Pro396Leu", "rs149047410"),
    ("JAK3", "Leu421Arg", "rs535740127"),
    ("JAK3", "Arg431Gln", "rs144953325"),
    ("JAK3", "Thr435Ile", "rs199706172"),
    ("JAK3", "Arg451Gln", "rs145751599"),
    ("JAK3", "Gly491Arg", "rs200112185"),
    ("JAK3", "Gln507Arg", "rs140690573"),
    ("JAK3", "Leu521Val", "rs55666418"),
    ("JAK3", "His529Arg", "rs142805245"),
    ("JAK3", "Ala572Val", "rs121913504"),
    ("JAK3", "Ala572Val", "rs121913504"),
    ("JAK3", "Ala572Val", "rs121913504"),
    ("JAK3", "Ala572Val", "rs121913504"),
    ("JAK3", "Ala572Val", "rs121913504"),
    ("JAK3", "Ala573Val", "rs2147686240"),
    ("JAK3", "Ala573Val", "rs2147686240"),
    ("JAK3", "Ala573Val", "rs2147686240"),
    ("JAK3", "Ala573Val", "rs2147686240"),
    ("JAK3", "Ala573Val", "rs2147686240"),
    ("JAK3", "Arg657Gln", "rs758959409"),
    ("JAK3", "Arg657Gln", "rs758959409"),
    ("JAK3", "Val722Ile", "rs3213409"),
    ("JAK3", "Val722Ile", "rs3213409"),
    ("JAK3", "Val722Ile", "rs3213409"),
    ("JAK3", "Val722Ile", "rs3213409"),
    ("JAK3", "Val722Ile", "rs3213409"),
    ("JAK3", "Val722Ile", "rs3213409"),
    ("JAK3", "Val722Ile", "rs3213409"),
    ("JAK3", "Arg799Cys", "rs201241352"),
    ("JAK3", "Gln827Glu", "rs144683649"),
    ("JAK3", "Ser835Cys", "rs201966394"),
    ("JAK3", "Arg840Cys", "rs200077579"),
    ("JAK3", "Ala877Val", "rs201869359"),
    ("JAK3", "His879Arg", "rs3179893"),
    ("JAK3", "Arg887Cys", "rs759015510"),
    ("JAK3", "Arg887His", "rs148688786"),
    ("JAK3", "Arg916Trp", "rs375807308"),
    ("JAK3", "Ala919Val", "rs767424476"),
    ("JAK3", "Arg925Cys", "rs149452625"),
    ("JAK3", "Arg925Ser", "rs149452625"),
    ("JAK3", "Ile1003Val", "rs137901277"),
    ("JAK3", "Tyr1023His", "rs145260622"),
    ("JAK3", "Leu1073Phe", "rs200580168"),
    ("JAK3", "Ala1090Thr", "rs144968714"),
    ("JAK3", "Glu1106Gly", "rs374152339"),
    ("TYK2", "Arg4His", "rs12720343"),
    ("TYK2", "Arg4Cys", "rs368801193"),
    ("TYK2", "Val15Ala", "rs144960992"),
    ("TYK2", "Val15Ile", "rs374780145"),
    ("TYK2", "Ala53Thr", "rs55762744"),
    ("TYK2", "Pro86Thr", "rs141466711"),
    ("TYK2", "Arg118Trp", "rs138742402"),
    ("TYK2", "Arg118Gln", "rs369530676"),
    ("TYK2", "Pro146Leu", "rs137904701"),
    ("TYK2", "Glu150Gly", "rs200718118"),
    ("TYK2", "Val164Met", "rs531355933"),
    ("TYK2", "Val164Leu", "rs531355933"),
    ("TYK2", "Ala165Thr", "rs370420133"),
    ("TYK2", "Leu194Phe", "rs373574074"),
    ("TYK2", "Arg197Cys", "rs376427383"),
    ("TYK2", "Arg197His", "rs12720263"),
    ("TYK2", "Arg217Cys", "rs754414161"),
    ("TYK2", "Arg231Trp", "rs201917359"),
    ("TYK2", "Arg231Trp", "rs201917359"),
    ("TYK2", "Arg235Trp", "rs767776173"),
    ("TYK2", "Arg243Trp", "rs200003143"),
    ("TYK2", "Arg265Gln", "rs775315394"),
    ("TYK2", "Arg269Cys", "rs368109068"),
    ("TYK2", "Arg274His", "rs371937276"),
    ("TYK2", "Arg294Trp", "rs367656874"),
    ("TYK2", "Arg294Gln", "rs201811978"),
    ("TYK2", "Val362Phe", "rs2304256"),
    ("TYK2", "Val362Phe", "rs2304256"),
    ("TYK2", "Val362Phe", "rs2304256"),
    ("TYK2", "Gly363Ser", "rs2304255"),
    ("TYK2", "Pro365Leu", "rs373929549"),
    ("TYK2", "Arg381Trp", "rs201240289"),
    ("TYK2", "Val386Met", "rs55956017"),
    ("TYK2", "Arg442Gln", "rs2304254"),
    ("TYK2", "Arg465Gln", "rs771922681"),
    ("TYK2", "Arg490His", "rs369615833"),
    ("TYK2", "Gly496Ser", "rs202072909"),
    ("TYK2", "Arg501Trp", "rs756944840"),
    ("TYK2", "Arg501Gln", "rs200643906"),
    ("TYK2", "Arg503Gln", "rs541104349"),
    ("TYK2", "Gly520Asp", "rs142576987"),
    ("TYK2", "Asp543Asn", "rs747616008"),
    ("TYK2", "Arg549His", "rs536181491"),
    ("TYK2", "Met564Arg", "rs72561470"),
    ("TYK2", "Pro571Ser", "rs867249402"),
    ("TYK2", "Val603Met", "rs140594440"),
    ("TYK2", "Gly634Glu", "rs189471343"),
    ("TYK2", "Ala652Thr", "rs144290894"),
    ("TYK2", "Thr668Met", "rs146430107"),
    ("TYK2", "Ile684Ser", "rs12720356"),
    ("TYK2", "Ile684Ser", "rs12720356"),
    ("TYK2", "Ile684Ser", "rs12720356"),
    ("TYK2", "Ile684Ser", "rs12720356"),
    ("TYK2", "Ile684Ser", "rs12720356"),
    ("TYK2", "Ile684Ser", "rs12720356"),
    ("TYK2", "Arg701Thr", "rs200791116"),
    ("TYK2", "Arg703Trp", "rs55882956"),
    ("TYK2", "Arg744Trp", "rs142676100"),
    ("TYK2", "Gly761Val", "rs201335603"),
    ("TYK2", "Gly761Val", "rs201335603"),
    ("TYK2", "Arg772Gln", "rs148823525"),
    ("TYK2", "Arg772Trp", "rs755089851"),
    ("TYK2", "Leu784Val", "rs759014184"),
    ("TYK2", "Asp810Val", "rs371070470"),
    ("TYK2", "Pro814Leu", "rs143743593"),
    ("TYK2", "Arg818His", "rs145437969"),
    ("TYK2", "Pro820His", "rs34046749"),
    ("TYK2", "Ala928Val", "rs35018800"),
    ("TYK2", "Ala928Val", "rs35018800"),
    ("TYK2", "Ala928Val", "rs35018800"),
    ("TYK2", "His993Tyr", "rs201397594"),
    ("TYK2", "Ala1004Gly", "rs202110875"),
    ("TYK2", "Pro1104Ala", "rs34536443"),
    ("TYK2", "Pro1104Ala", "rs34536443"),
    ("TYK2", "Pro1104Ala", "rs34536443"),
    ("TYK2", "Glu1163Gly", "rs55886939"),
    ("TYK2", "Cys1187Tyr", "rs200932305")
    ]


# Function to query ClinVar for each rsID
def search_clinvar(rsid):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "clinvar", "term": rsid, "retmode": "json"}

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        result = response.json()
        return result.get("esearchresult", {}).get("idlist", [])
    except requests.RequestException as e:
        print(f"Error querying {rsid}: {e}")
        return []


# Function to fetch clinical significance from ClinVar
def fetch_clinvar_variant(variant_id):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": "clinvar", "id": variant_id, "retmode": "json"}

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        result = response.json()

        for key in result["result"]:
            if key != "uids" and "germline_classification" in result["result"][key]:
                classification = result["result"][key]["germline_classification"]
                if "description" in classification:
                    return classification["description"]  # Extract correct significance

    except requests.RequestException as e:
        print(f"Error fetching ClinVar data for {variant_id}: {e}")

    return "Unknown"


# Store all variant results
variant_significance = []

# Iterate through rsIDs and query ClinVar
for gene, amino_acid_change, rsid in extracted_variants_list:
    print(f"Processing {rsid}...")

    if rsid == "N/A":
        significance = "Not Available"
    else:
        variant_ids = search_clinvar(rsid)

        if not variant_ids:
            significance = "No Match Found"
        else:
            significance = "Unknown"
            for vid in variant_ids:
                significance = fetch_clinvar_variant(vid)
                if significance != "Unkown":
                    print("FOUND ONE")
                if significance and significance != "Unknown":
                    break  # Stop if we get a valid significance

    variant_significance.append((gene, amino_acid_change, rsid, significance))
    time.sleep(1)  # Prevent API rate limiting

# Convert to DataFrame
df_variants = pd.DataFrame(variant_significance,
                           columns=["Gene", "Amino Acid Substitution", "rsID", "Clinical Significance"])

# Save full list of variants with ClinVar results
df_variants.to_csv(full_output_path, sep="\t", index=False)
print(f"Full results saved to {full_output_path}")

# Filter for benign variants only
df_benign = df_variants[df_variants["Clinical Significance"].str.contains("benign", case=False, na=False)]

# Save benign-only variants
df_benign.to_csv(benign_output_path, sep="\t", index=False)
print(f"Benign results saved to {benign_output_path}")