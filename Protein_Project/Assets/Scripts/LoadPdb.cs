using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System;

public class LoadPdb : MonoBehaviour
{
    public class CAtom
    {
        public int SerialNumber;
        public int ResID;
        public string Name;
        public string ResName;
        public float Xcoordinate;
        public float Ycoordinate;
        public float Zcoordinate;
        public float BFactor;
        public GameObject gObj;
        public List<CAtom> bonds = new List<CAtom>();
    }
    List<CAtom> LAtoms = new List<CAtom>();
    public GameObject Sph;
    void LoadFromPdb()
    {
        string m_Path = Application.dataPath;
        string strPdb = "1afb.pdb";
        StreamReader sr = new StreamReader(m_Path + "/" + strPdb);
        string str;

        while (true)
        {

            CAtom atom = new CAtom();
            str = sr.ReadLine();
            if (str == null)
            {
                break;
            }
            if (str.StartsWith("ATOM"))
            {

                atom.SerialNumber = int.Parse(str.Substring(7, 5).Trim());
                atom.Name = str.Substring(13, 4).Trim();
                atom.ResName = str.Substring(17, 3).Trim();
                atom.ResID = int.Parse(str.Substring(23, 4).Trim());
                atom.Xcoordinate = float.Parse(str.Substring(31, 8).Trim());
                atom.Ycoordinate = float.Parse(str.Substring(39, 8).Trim());
                atom.Zcoordinate = float.Parse(str.Substring(47, 8).Trim());
                atom.BFactor = float.Parse(str.Substring(61, 5).Trim());
                Vector3 pos = new Vector3(atom.Xcoordinate, atom.Ycoordinate, atom.Zcoordinate);
                atom.gObj = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                atom.gObj.transform.position = pos;
                LAtoms.Add(atom);

            }
        }
        sr.Close();
    }
    void PrintFile()
    {
        foreach (var item in LAtoms)
        {
            print("ATOM" + "      " + item.SerialNumber + "    " + item.Name + "    " + item.ResName + "    " + item.ResID + "    " + item.Xcoordinate + "    " + item.Ycoordinate + "    " + item.Zcoordinate + "     " + item.BFactor);
        }
    }
    void ColorAtoms()
    {
        for (int i = 0; i < LAtoms.Count; i++)
        {
            Renderer Rnd = LAtoms[i].gObj.GetComponent<Renderer>();
            if (LAtoms[i].Name == "O")
            {
                Rnd.material.SetColor("_Color", Color.red);
            }
            if (LAtoms[i].Name == "C")
            {
                Rnd.material.SetColor("_Color", Color.black);
            }
            if (LAtoms[i].Name == "N")
            {
                Rnd.material.SetColor("_Color", Color.blue);
            }
            if (LAtoms[i].Name == "H")
            {
                Rnd.material.SetColor("_Color", Color.gray);
            }

        }
    }
    void HowManyAtoms()
    {
        int countN = 0;
        int countH = 0;
        int countO = 0;
        int countC = 0;

        for (int i = 0; i < LAtoms.Count; i++)
        {
            if (LAtoms[i].Name == "N")
            {
                countN++;
            }
            if (LAtoms[i].Name == "C")
            {
                countC++;
            }
            if (LAtoms[i].Name == "O")
            {
                countO++;
            }
            if (LAtoms[i].Name == "H")
            {
                countH++;
            }
        }
        print("N ATOMS = " + countN);
        print("C ATOMS = " + countC);
        print("H ATOMS = " + countH);
        print("O ATOMS = " + countO);
    }
    void AppendBond()
    {
        for (int i = 0; i < LAtoms.Count; i++)
        {
            for (int j = i + 1; j < LAtoms.Count; j++)
            {
                float dx = LAtoms[i].gObj.transform.position.x - LAtoms[j].gObj.transform.position.x;
                float dz = LAtoms[i].gObj.transform.position.z - LAtoms[j].gObj.transform.position.z;
                float dy = LAtoms[i].gObj.transform.position.y - LAtoms[j].gObj.transform.position.y;

                double dis = Mathf.Sqrt((dx * dx) + (dy * dy) + (dz * dz));

                if (dis < 1.6)
                {
                    LAtoms[i].bonds.Add(LAtoms[j]);
                    LAtoms[j].bonds.Add(LAtoms[i]);
                }


            }
        }

    }
    void AvgBondLength() // calcluate the average bond length of atoms "N" and the ResID = 73 
    {
        for (int i = 0; i < LAtoms.Count; i++)
        {
            if (LAtoms[i].Name == "N" && LAtoms[i].ResID == 73)
            {
                float sum = 0;
                for (int j = 0; j < LAtoms[i].bonds.Count; j++)
                {
                    CAtom BondenAtom = LAtoms[i].bonds[j];
                    float dx = BondenAtom.gObj.transform.position.x - LAtoms[i].gObj.transform.position.x;
                    float dz = BondenAtom.gObj.transform.position.z - LAtoms[i].gObj.transform.position.z;
                    float dy = BondenAtom.gObj.transform.position.y - LAtoms[i].gObj.transform.position.y;

                    float dist = Mathf.Sqrt(dx * dx + dy * dy + dz * dz);
                    sum += dist;
                }

                float avg = sum / LAtoms[i].bonds.Count;
                print("Average bond length for  the atom ('N') of ResID('73')   = " + avg);

            }
        }
    }
    void AvgBondSpherical(float Xcenter, float Ycenter, float Zcenter, float R) //calculate the avg bond length of atoms that located in the given circle
    {
        for (int i = 0; i < LAtoms.Count; i++)
        {
            float CalcX = (LAtoms[i].gObj.transform.position.x - Xcenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);
            float CalcY = (LAtoms[i].gObj.transform.position.y - Ycenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);
            float CalcZ = (LAtoms[i].gObj.transform.position.z - Zcenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);

            float CalcDistance = Mathf.Sqrt(CalcX + CalcY + CalcZ); 

            if (CalcDistance < R)
            {

                float sum = 0;
                for (int j = 0; j < LAtoms[i].bonds.Count; j++)
                {
                    CAtom BondenAtom = LAtoms[i].bonds[j];
                    float dx = BondenAtom.gObj.transform.position.x - LAtoms[i].gObj.transform.position.x;
                    float dz = BondenAtom.gObj.transform.position.z - LAtoms[i].gObj.transform.position.z;
                    float dy = BondenAtom.gObj.transform.position.y - LAtoms[i].gObj.transform.position.y;

                    float dist = Mathf.Sqrt(dx * dx + dy * dy + dz * dz);
                    sum += dist;
                }

                float avg = sum / LAtoms[i].bonds.Count;
                print("Average bond length for  the atom ( " + LAtoms[i].Name + " ) " + " =  " + avg);
            }


        }
    }
    void NumberOfAtomsSpherical(float Xcenter, float Ycenter, float Zcenter, float R) // calculate the number of atoms that located in a given circle
    {
        int CounterAtoms = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            float CalcX = (LAtoms[i].gObj.transform.position.x - Xcenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);
            float CalcY = (LAtoms[i].gObj.transform.position.y - Ycenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);
            float CalcZ = (LAtoms[i].gObj.transform.position.z - Zcenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);

            float CalcDistance = Mathf.Sqrt(CalcX + CalcY + CalcZ);

            if ( CalcDistance < R )
            {
                CounterAtoms += 1;
            }
        
        }
        print("Total Number of Atoms inside the Spherical = " + CounterAtoms);
    }
    void NumberOfResiduesSpherical(float Xcenter, float Ycenter, float Zcenter, float R) // calculate the number of Residues that located in a given circle
    {
        int CounterRes = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            if (LAtoms[i].Name == "CA" )
            {
                float CalcX = (LAtoms[i].gObj.transform.position.x - Xcenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);
                float CalcY = (LAtoms[i].gObj.transform.position.y - Ycenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);
                float CalcZ = (LAtoms[i].gObj.transform.position.z - Zcenter) * (LAtoms[i].gObj.transform.position.x - Xcenter);

                float CalcDistance = Mathf.Sqrt(CalcX + CalcY + CalcZ);

                if (CalcDistance < R)
                {
                    CounterRes += 1;
                }
            }

        }
        print("Total Number of Residues inside the Spherical = " + CounterRes); // calculate the total bond lengths for the Residues (73:83)
    }
    void TotalBonds()
    {
        float sum = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            if (LAtoms[i].Name == "CA" && LAtoms[i].ResID >= 73 && LAtoms[i].ResID <= 82)
            {
                for (int j = i + 1; j < LAtoms.Count; j++)
                {
                    if (LAtoms[j].Name == "CA" && LAtoms[j].ResID == (LAtoms[i].ResID + 1))
                    {
                        CAtom BondenAtom = LAtoms[j];
                        float dx = BondenAtom.gObj.transform.position.x - LAtoms[i].gObj.transform.position.x;
                        float dy = BondenAtom.gObj.transform.position.y - LAtoms[i].gObj.transform.position.y;
                        float dz = BondenAtom.gObj.transform.position.z - LAtoms[i].gObj.transform.position.z;

                        float dist = Mathf.Sqrt(dx * dx + dy * dy + dz * dz);
                        sum += dist;
                    }
                }
            }
        }
        print("Total bond lengths for the Residues (73:83)  = " + sum);
    }
    void AvgBfactorAll() //calculate the avg Bfactor of all atoms of pdb file
    {
        float sum = 0;
        float avg = 0;
        for(int i=0; i< LAtoms.Count ;i++)
        {
            sum += LAtoms[i].BFactor;
        }
        avg = sum / LAtoms.Count;
        print("avg Bfactor of all atoms   = " + avg);
    }
    void AvgBfactorSpecific() //calculate the avg Bfactor of  the residues (35 : 60)
    {
        float sum = 0;
        float avg = 0;
        int counter = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            if (LAtoms[i].ResID >= 35 && LAtoms[i].ResID <= 60)
            {

                counter += 1;
                sum += LAtoms[i].BFactor;
            }
        }
        avg = sum / counter;
        print("avg  Bfactor of  the residues (35 : 60)    = " + avg);
    }
    void AvgBfactorRegion(float xMin, float xMax,float yMin, float yMax,float zMin, float zMax) //calculate the avg Bfactor of all atoms at a given region
    {
        float sum = 0;
        int count = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            Vector3 pos = LAtoms[i].gObj.transform.position;
            if (pos.x >= xMin && pos.x <= xMax && pos.y >= yMin && pos.y <= yMax && pos.z >= zMin && pos.z <= zMax)
            {
                sum += LAtoms[i].BFactor;
                count++;
            }
        }
        float avg = sum / count;
        print("Average Bfactor of region  ( Xmin = " + xMin +" , " + "Xmax = "+xMax + " , " +"Ymin = "+yMax + " , " +"Ymax = "+yMax+" , "+"Zmin = "+zMin+" , "+"Zmax = "+zMax+" )"+" = "+ avg);

    }
    void AvgBfactorRegionCA(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax) //calculate the average Bfactor of region CA atom
    {
        float sum = 0;
        int count = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            Vector3 pos = LAtoms[i].gObj.transform.position;
            if (LAtoms[i].Name=="CA" &&  pos.x >= xMin && pos.x <= xMax && pos.y >= yMin && pos.y <= yMax && pos.z >= zMin && pos.z <= zMax)
            {
                sum += LAtoms[i].BFactor;
                count++;
            }
        }

        float avg = sum / count;
        print("Average Bfactor of region CA ATOM  ( Xmin = " + xMin + " , " + "Xmax = " + xMax + " , " + "Ymin = " + yMax + " , " + "Ymax = " + yMax + " , " + "Zmin = " + zMin + " , " + "Zmax = " + zMax + " )" + " = " + avg);

    }
    void DrawBond() //Draw bonds between all atoms
    {
        for (int i = 0; i < LAtoms.Count; i++)
        {
            for (int j = 0; j < LAtoms[i].bonds.Count; j++)
            {
                Debug.DrawLine(LAtoms[i].gObj.transform.position, LAtoms[i].bonds[j].gObj.transform.position, Color.black, 2.5f);
            }
        }
    }
    void RotationX() //Rotate all atoms about x-axis at theta = 0.1
    {
        float nX, nY, nZ;
        float th = 0.1f;
        for(int i=0;i<LAtoms.Count;i++)
        {
            Vector3 pos = LAtoms[i].gObj.transform.position;
            nY = pos.y * Mathf.Cos(th) - pos.z * Mathf.Sin(th);
            nZ = pos.y * Mathf.Sin(th) + pos.z * Mathf.Cos(th);
            nX = pos.x;
            LAtoms[i].gObj.transform.position = new Vector3(nX, nY, nZ);
        }
    }
    void RotationZ() //Rotate all atoms about z-axis at theta = 0.1
    {
        float nX, nY, nZ;
        float th = 0.1f;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            Vector3 pos = LAtoms[i].gObj.transform.position;
            nY = pos.x * Mathf.Sin(th) + pos.y * Mathf.Cos(th);
            nX = pos.x * Mathf.Cos(th) - pos.y * Mathf.Sin(th);
            nZ = pos.z;
            LAtoms[i].gObj.transform.position = new Vector3(nX, nY, nZ);
        }
    }
    void RotationY() //Rotate all atoms about y-axis at theta = 0.1
    {
        float nX, nY, nZ;
        float th = 0.1f;
        for (int i = 0; i < LAtoms.Count; i++)
        {
            Vector3 pos = LAtoms[i].gObj.transform.position;
            nY = pos.y;
            nX = pos.z * Mathf.Sin(th) + pos.x * Mathf.Cos(th);  
            nZ = pos.z * Mathf.Cos(th) - pos.x * Mathf.Sin(th);
            LAtoms[i].gObj.transform.position = new Vector3(nX, nY, nZ);
        }
    }
    void RotateArbitaty(Vector3 p1,Vector3 p2) //Rotate all atoms arbitatry at theta = 0.1
    {
        float nX, nY, nZ;
        float th = 0.1f;
        float dx = p2.x - p1.x;
        float dy = p2.y - p1.y;
        float dz = p2.z - p1.z;

        float theta = (float)Mathf.Atan2(dy, dx);
        float phi = (float)Mathf.Atan2(Mathf.Sqrt(dx * dx + dy * dy), dz);

        for(int i=0;i<LAtoms.Count;i++)
        {
            Vector3 pos = LAtoms[i].gObj.transform.position;
            pos.x -= p1.x;
            pos.y -= p1.y;
            pos.z -= p1.z;

            // Rz^-1(theta)
            nY = pos.y * Mathf.Sin(-theta) + pos.x * Mathf.Cos(-theta);
            nX = pos.x * Mathf.Cos(-theta) - pos.y * Mathf.Sin(-theta);
            nZ = pos.z;
            pos = new Vector3(nX, nY, nZ);

            // Ry^-1(phi)
            nZ = pos.z * Mathf.Cos(-phi) - pos.x * Mathf.Sin(-phi);
            nX= pos.z * Mathf.Sin(-phi) + pos.x * Mathf.Cos(-phi);
            nY = pos.y;
            pos = new Vector3(nX, nY, nZ);


            // Rz(th)
            nX = pos.x * Mathf.Cos(th) - pos.y * Mathf.Sin(th);
            nY= pos.y * Mathf.Sin(th) + pos.x * Mathf.Cos(th);
            nZ = pos.z;
            pos = new Vector3(nX, nY, nZ);

            // Ry(phi)
            nZ = pos.z * Mathf.Cos(phi) - pos.x * Mathf.Sin(phi);
            nX = pos.z * Mathf.Sin(phi) + pos.x * Mathf.Cos(phi);
            nY = pos.y;
            pos = new Vector3(nX, nY, nZ);

            //// Rz(theta)
            nY = pos.y * Mathf.Sin(theta) + pos.x * Mathf.Cos(theta);
            nX = pos.x * Mathf.Cos(theta) - pos.y * Mathf.Sin(theta);
            nZ = pos.z;
            pos = new Vector3(nX, nY, nZ);

            pos.x += p1.x;
            pos.y += p1.y;
            pos.z += p1.z;
            LAtoms[i].gObj.transform.position = pos;

        }

    }
    void AvgBFactorBetweenTwoBox(float xMin, float xMax,
                   float yMin, float yMax,
                   float zMin, float zMax,
                   float xMin2, float xMax2,
                   float yMin2, float yMax2,
                   float zMin2, float zMax2)  //calculate the average BFactor for atoms between two given Boxes 
    {
        float sum = 0;
        int counter = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {

            if (LAtoms[i].gObj.transform.position.x >= xMin && LAtoms[i].gObj.transform.position.x <= xMax
                    && LAtoms[i].gObj.transform.position.y >= yMin && LAtoms[i].gObj.transform.position.y <= yMax
                    && LAtoms[i].gObj.transform.position.z >= zMin && LAtoms[i].gObj.transform.position.z <= zMax
                    && LAtoms[i].gObj.transform.position.x >= xMin2 && LAtoms[i].gObj.transform.position.x <= xMax2
                    && LAtoms[i].gObj.transform.position.y >= yMin2 && LAtoms[i].gObj.transform.position.y <= yMax2
                    && LAtoms[i].gObj.transform.position.z >= zMin2 && LAtoms[i].gObj.transform.position.z <= zMax2
                )
            {
                sum += LAtoms[i].BFactor;
                counter++;
            }
        }
        float avg = sum / counter;
        print("the average BFactor for atoms between two given Box region  = " + avg);
    }
    void AtomsInsSpherical(float xCenter1, float yCenter1, float zCeneter1, float Raduis1,
        float xCenter2, float yCenter2, float zCeneter2, float Raduis2) //calculate the AVG bond length of all atoms that between two given speiral
    {
        float dx1;
        float dz1;
        float dy1;
        float dx2;
        float dz2;
        float dy2;
        float Counter = 0;
        float TotalBonds = 0;
        for (int i = 0; i < LAtoms.Count; i++)
        {


            dx1 = LAtoms[i].gObj.transform.position.x - xCenter1;
            dy1 = LAtoms[i].gObj.transform.position.y - yCenter1;
            dz1 = LAtoms[i].gObj.transform.position.z - zCeneter1;

            float dist = Mathf.Sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
            if (dist > Raduis1)
            {

                dx2 = LAtoms[i].gObj.transform.position.x - xCenter2;
                dy2 = LAtoms[i].gObj.transform.position.y - yCenter2;
                dz2 = LAtoms[i].gObj.transform.position.z - zCeneter2;

                float dist2 = Mathf.Sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                if (dist < Raduis2)
                {
                    float sum = 0;
                    for (int j = 0; j < LAtoms[i].bonds.Count; j++)
                    {
                        CAtom other = LAtoms[i].bonds[j];
                        float dx3 = other.gObj.transform.position.x - LAtoms[i].gObj.transform.position.x;
                        float dy3 = other.gObj.transform.position.y - LAtoms[i].gObj.transform.position.y;
                        float dz3 = other.gObj.transform.position.z - LAtoms[i].gObj.transform.position.z;

                        float dist3 = Mathf.Sqrt(dx3 * dx3 + dy3 * dy3 + dz3 * dz3);
                        sum += dist3;

                    }

                    float avg = sum / LAtoms[i].bonds.Count;
                    print("AVG Bond length  atoms between  the  two speiral is   = " + avg);
                }
            }
        }

    }   
    void Start()
    {
        LoadFromPdb();
        //PrintFile()
        //HowManyAtoms();
        AppendBond();
        //AvgBondLength(); -->Done
        //TotalBonds(); --> Done
        // AvgBfactorAll(); -->Done
        //AvgBfactorSpecific(); --> Done
        //AvgBfactorRegion(-14.114f, 1.630f, 20.631f, 28.122f, 31.379f, 41.703f); -->Done
        //AvgBfactorRegionCA(-14.114f, 1.630f, 20.631f, 28.122f, 31.379f, 41.703f); -->Done
        //AvgBondSpherical(-14.114f,28.122f,41.703f,0.1f); -->Done
        //NumberOfAtomsSpherical(17.554f, 21.744f, 22.807f, 0.1f); -->Done
        //NumberOfResiduesSpherical(17.554f, 21.744f, 22.807f, 0.1f); -->Done
        
    }

    // Update is called once per frame
    void Update()
    {
        //ColorAtoms();
        //RotationX();
        //RotationZ();
        //RotationY();
        Vector3 p1 = new Vector3(-10.952f, 24.905f, 44.032f);
        Vector3 p2 = new Vector3(36.834f, -8.661f , 24.348f);
        //RotateArbitaty(p1,p2);
        // DrawBond();
    }
}

