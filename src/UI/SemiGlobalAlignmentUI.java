package UI;

import semiglobalalignment.SemiGlobalAlignment;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Kiky
 */
//import localalignment.SemiGlobalAlignmentUI;

import javax.swing.JOptionPane;
public class SemiGlobalAlignmentUI extends javax.swing.JFrame {
    //public static SemiGlobalAlignmentUI main= new SemiGlobalAlignmentUI();
    /**
     * Creates new form SemiGlobalAlignmentUI
     */
    public SemiGlobalAlignmentUI() {
        initComponents();
	errorQ.setVisible(false);
	errorD.setVisible(false);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        sequenceQ = new java.awt.TextField();
        labelSeqQ = new java.awt.Label();
        sequenceD = new java.awt.TextField();
        labelSeqD = new java.awt.Label();
        Align = new javax.swing.JButton();
        properties = new java.awt.Button();
        errorQ = new javax.swing.JLabel();
        errorD = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Semi Global Alignment");
        setFocusTraversalPolicyProvider(true);
        setName("Local Alignment"); // NOI18N

        sequenceQ.setForeground(new java.awt.Color(0, 0, 0));
        sequenceQ.setText("sequence q");
        sequenceQ.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                sequenceQMouseClicked(evt);
            }
        });
        sequenceQ.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                sequenceQFocusLost(evt);
            }
        });
        sequenceQ.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sequenceQActionPerformed(evt);
            }
        });

        labelSeqQ.setFont(new java.awt.Font("Dialog", 2, 12)); // NOI18N
        labelSeqQ.setText("Please input protein sequence q (separate by a space)");

        sequenceD.setForeground(new java.awt.Color(0, 0, 0));
        sequenceD.setText("sequence d");
        sequenceD.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                sequenceDMouseClicked(evt);
            }
        });
        sequenceD.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                sequenceDFocusLost(evt);
            }
        });

        labelSeqD.setFont(new java.awt.Font("Dialog", 2, 12)); // NOI18N
        labelSeqD.setText("Please input protein sequence d  (separate by a space)");

        Align.setFont(new java.awt.Font("Tahoma", 1, 11)); // NOI18N
        Align.setLabel("Align");
        Align.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                AlignMouseClicked(evt);
            }
        });

        properties.setLabel("Properties");
        properties.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                propertiesActionPerformed(evt);
            }
        });

        errorQ.setForeground(new java.awt.Color(204, 0, 0));
        errorQ.setText("this field can not be empty");

        errorD.setForeground(new java.awt.Color(204, 0, 0));
        errorD.setText("this field can not be empty");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(errorD)
                        .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addComponent(labelSeqD, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(sequenceD, javax.swing.GroupLayout.PREFERRED_SIZE, 264, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(errorQ)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addComponent(labelSeqQ, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(sequenceQ, javax.swing.GroupLayout.PREFERRED_SIZE, 261, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addContainerGap(48, Short.MAX_VALUE))))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(properties, javax.swing.GroupLayout.PREFERRED_SIZE, 92, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(22, 22, 22)
                .addComponent(Align, javax.swing.GroupLayout.PREFERRED_SIZE, 69, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(39, 39, 39)
                .addComponent(labelSeqQ, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(14, 14, 14)
                .addComponent(sequenceQ, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(errorQ)
                .addGap(25, 25, 25)
                .addComponent(labelSeqD, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(sequenceD, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(errorD)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 33, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(Align, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(properties, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );

        setBounds(200, 100, 379, 332);
    }// </editor-fold>//GEN-END:initComponents

    private void sequenceQActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sequenceQActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_sequenceQActionPerformed

    private void sequenceDMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_sequenceDMouseClicked
        if (sequenceD.isFocusOwner()==true && sequenceD.getText().equals("sequence d"))
		sequenceD.setText(null);
	errorD.setVisible(false);// TODO add your handling code here:
    }//GEN-LAST:event_sequenceDMouseClicked

    private void sequenceQMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_sequenceQMouseClicked
        if (sequenceQ.isFocusOwner()==true && sequenceQ.getText().equals("sequence q"))
		sequenceQ.setText(null);
	errorQ.setVisible(false);// TODO add your handling code here:
    }//GEN-LAST:event_sequenceQMouseClicked

    private void sequenceQFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_sequenceQFocusLost
        if (sequenceQ.getText().equals(""))
	    sequenceQ.setText("sequence q");// TODO add your handling code here:
    }//GEN-LAST:event_sequenceQFocusLost

    private void sequenceDFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_sequenceDFocusLost
        if (sequenceD.getText().equals(""))
	    sequenceD.setText("sequence d");// TODO add your handling code here:
    }//GEN-LAST:event_sequenceDFocusLost

    private void propertiesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_propertiesActionPerformed
         this.setEnabled(false);
	 Properties prop =new Properties();
	 prop.setVisible(true);
	// TODO add your handling code here:
    }//GEN-LAST:event_propertiesActionPerformed

    private void AlignMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_AlignMouseClicked
        if(sequenceQ.getText().equals("sequence q"))
	    errorQ.setVisible(true);
	if(sequenceD.getText().equals("sequence d"))
	    errorD.setVisible(true);
	if (errorQ.isVisible()==false && errorD.isVisible()==false){
	    semiglobalalignment.SemiGlobalAlignment.fungsi.dataD=sequenceD.getText().toLowerCase();
	    semiglobalalignment.SemiGlobalAlignment.fungsi.dataQ=sequenceQ.getText().toLowerCase();
	    if (SemiGlobalAlignment.fungsi.pilihangap==1){
                SemiGlobalAlignment.fungsi.matrikstest();
                if (SemiGlobalAlignment.fungsi.pilihanmatriks==1){
                    SemiGlobalAlignment.fungsi.algoritmaPAM250();
                    SemiGlobalAlignment.fungsi.alignment();
                }
                else if (SemiGlobalAlignment.fungsi.pilihanmatriks==2){
                    SemiGlobalAlignment.fungsi.algoritmaBlosum62();
                    SemiGlobalAlignment.fungsi.alignment();
                }
		//else
		//    JOptionPane.showMessageDialog(null,"matrix has not been initialized!");
            }
            else if (SemiGlobalAlignment.fungsi.pilihangap==2){
                SemiGlobalAlignment.fungsi.matrikstest();
                if (SemiGlobalAlignment.fungsi.pilihanmatriks==1){
                    SemiGlobalAlignment.fungsi.algoritmaPAM250();
                    SemiGlobalAlignment.fungsi.alignment();
                }
                else if (SemiGlobalAlignment.fungsi.pilihanmatriks==2){
                    SemiGlobalAlignment.fungsi.algoritmaBlosum62();
                    SemiGlobalAlignment.fungsi.alignment();
                }
            }
	    else
		JOptionPane.showMessageDialog(null,"gap and matrix has not been initialized!");
	}// TODO add your handling code here:
    }//GEN-LAST:event_AlignMouseClicked

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(SemiGlobalAlignmentUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(SemiGlobalAlignmentUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(SemiGlobalAlignmentUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(SemiGlobalAlignmentUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
		//main.setVisible(true);
	    }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton Align;
    private javax.swing.JLabel errorD;
    private javax.swing.JLabel errorQ;
    private java.awt.Label labelSeqD;
    private java.awt.Label labelSeqQ;
    private java.awt.Button properties;
    private java.awt.TextField sequenceD;
    private java.awt.TextField sequenceQ;
    // End of variables declaration//GEN-END:variables
}
