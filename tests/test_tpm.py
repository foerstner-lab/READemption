from reademptionlib.genewisequanti import GeneWiseOverview
import pandas as pd


def test_tpm_calculation():
    raw = pd.DataFrame(
        data={
            "Orientation of counted reads relative to the strand location of the annotation": [
                "sense",
                "sense",
                "sense",
            ],
            "Sequence name": ["NC_016810.1", "NC_016810.1", "NC_016810.1"],
            "Source": ["RefSeq", "RefSeq", "RefSeq"],
            "Feature": ["CDS", "CDS", "CDS"],
            "Start": ["1", "51", "101"],
            "End": ["25", "100", "200"],
            "Score": [".", ".", "."],
            "Strand": ["+", "+", "+"],
            "Frame": ["0", "0", "0"],
            "Attributes": [
                "ID=cds0;Name=YP_005179941.1",
                "ID=cds1;Name=YP_005179942.1",
                "ID=cds2;Name=YP_005179943.1",
            ],
            "InSPI2_R1": ["6", "10", "80"],
            "InSPI2_R2": ["10", "20", "20"],
        }
    )
    expected_tpm = pd.DataFrame(
        data={
            "Orientation of counted reads relative to the strand location of the annotation": [
                "sense",
                "sense",
                "sense",
            ],
            "Sequence name": ["NC_016810.1", "NC_016810.1", "NC_016810.1"],
            "Source": ["RefSeq", "RefSeq", "RefSeq"],
            "Feature": ["CDS", "CDS", "CDS"],
            "Start": ["1", "51", "101"],
            "End": ["25", "100", "200"],
            "Score": [".", ".", "."],
            "Strand": ["+", "+", "+"],
            "Frame": ["0", "0", "0"],
            "Attributes": [
                "ID=cds0;Name=YP_005179941.1",
                "ID=cds1;Name=YP_005179942.1",
                "ID=cds2;Name=YP_005179943.1",
            ],
            "InSPI2_R1": ["193548.387097", "161290.322581", "645161.290323"],
            "InSPI2_R2": ["400000.0", "400000.0", "200000.0"],
        }
    )
    expected_tpm["InSPI2_R1"] = expected_tpm["InSPI2_R1"].astype("float64")
    expected_tpm["InSPI2_R2"] = expected_tpm["InSPI2_R2"].astype("float64")
    calculated_tpm = GeneWiseOverview._calculate_tpm(
        "self", raw, ["InSPI2_R1", "InSPI2_R2"]
    )
    pd.testing.assert_frame_equal(expected_tpm, calculated_tpm)
